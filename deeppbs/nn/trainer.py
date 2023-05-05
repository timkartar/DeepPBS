# builtin modules
import logging
import json
from os.path import join as ospj
from collections import OrderedDict
import copy, sys

# third party modules
import torch
from torch.optim.lr_scheduler import ReduceLROnPlateau, OneCycleLR, ExponentialLR
from torch_geometric.nn import DataParallel
import torch.nn.functional as F

# deeppbs modules
from deeppbs.nn.utils import classWeights
from deeppbs.nn import processBatch
from deeppbs.nn.metrics import reportMetrics

class Scheduler(object):
    def __init__(self, scheduler):
        self.epoch = 0
        self.scheduler = scheduler
        self.history = {
            "loss": 0,
            "batch_count": 0
        }
    
    def step(self, epoch, loss, **kwargs):
        # per-epoch schedulers
        new_epoch = (epoch > self.epoch)
        if new_epoch:
            # we are in a new epoch, update per-epoch schedulers
            self.epoch = epoch
            if isinstance(self.scheduler, ReduceLROnPlateau):
                mean_loss = self.history['loss']/self.history["batch_count"]
                self.scheduler.step(mean_loss)
                self.history["loss"] = 0
                self.history["batch_count"] = 0
            elif isinstance(self.scheduler, ExponentialLR):
                self.scheduler.step()
        
        # per-batch schedulers
        if isinstance(self.scheduler, OneCycleLR):
            self.scheduler.step()
        elif isinstance(self.scheduler, ReduceLROnPlateau):
            self.history["loss"] += loss
            self.history["batch_count"] += 1

class Trainer(object):
    def __init__(self, model, nc, optimizer, criterion, 
            device='cpu', scheduler=None, evaluator=None, 
            writer=None, checkpoint_path='.', quiet=True,
            ic_loss_weight=1, mse_loss_weight=1, ce_loss_weight=1, kl_loss_weight=1
            ):
        # parameters
        self.model = model
        self.nc = nc
        self.optimizer = optimizer
        self.criterion = criterion
        self.kldiv = torch.nn.KLDivLoss(reduction = 'none')#, log_target=True) #Information content same as kldiv relative to 
        #self.kldivloss = torch.nn.KLDivLoss(reduction = 'none')
        self.ce_loss_weight = ce_loss_weight
        self.ic_loss_weight = ic_loss_weight
        self.mse_loss_weight = mse_loss_weight
        self.kl_loss_weight = kl_loss_weight
        self.evaluator = evaluator
        self.writer = writer
        self.device = device
        self.quiet = quiet
        self.checkpoint_path = checkpoint_path
        
        # variables to track training progress
        self.best_state = None
        self.best_state_metric = None
        self.best_epoch = None
        
        # set up scheduler
        if scheduler is not None:
            scheduler = Scheduler(scheduler)
        self.scheduler = scheduler
        
        # get model name
        if isinstance(self.model, DataParallel):
            self.model_name = self.model.module.name
        else:
            self.model_name = self.model.name
        
        # history
        self.metrics_history = {'epochs': []}
    
   
    def train(self, nepochs, dataset,
        validation_dataset=None, batch_loss_every=4, eval_every=2, debug=False,
        checkpoint_every=None, optimizer_kwargs={}, scheduler_kwargs={},
        best_state_metric=None, best_state_metric_threshold=None, 
        best_state_metric_dataset='validation', best_state_metric_goal='max', params_to_write=None,
        metrics_calculation="average_batches", use_mask=True
    ):
        # begin training
        if not self.quiet:
            logging.info("Beginning Training ({} epochs)".format(nepochs))
        
        if debug:
            mem_stats = {
                "current": [],
                "peak": [],
                "epoch_start": []
            }
        
        if best_state_metric_goal == 'max':
            self.best_state_metric = -999999
        else:
            self.best_state_metric = 999999
        
        batch_count = 0
        first_epoch = True
        for epoch in range(nepochs):
            # set model to training mode
            self.model.train()
            
            # forward + backward + update
            epoch_loss = 0
            n = 0
            for batch in dataset:
                oom = False
                # update the model weights
                batch_data = processBatch(self.device, batch)
                #batch, y, mask = batch_data['batch'], batch_data['y'], batch_data['mask']
                
                # check for OOM errors
                #loss = self.optimizer_step(batch_data, **optimizer_kwargs)
                #try:
                loss = self.optimizer_step(batch_data, epoch=epoch, **optimizer_kwargs)
                '''
                except RuntimeError as e: # out of memory
                    logging.info("Runtime error -- skipping batch.")
                    logging.debug("Error at loss computation.", exc_info=e)
                    oom = True
                '''
                if oom:
                    continue
                    
                # update scheduler
                if self.scheduler is not None:
                    self.scheduler.step(epoch, loss, **scheduler_kwargs)
                
                # write batch-level stats
                if batch_count % batch_loss_every  == 0:
                    if self.writer:
                        self.writer.add_scalar("train/batch_loss", loss, batch_count)
                
                ######### adding parameters of model to tensorboard ######
                if (params_to_write is not None) and self.writer:
                    for name, param in self.model.named_parameters():
                        if (name in params_to_write) and param.requires_grad:
                            if(param.data.cpu().numpy().flatten().shape[0] == 1):
                                self.writer.add_scalar(name, param.data.cpu().numpy()[0], batch_count)
                                self.writer.add_scalar(name + "_grad", param.grad.cpu().numpy()[0], batch_count)
                            elif(param.data.cpu().numpy().flatten().shape[0] <= 4):
                                for i in range(1,param.data.cpu().numpy().flatten().shape[0] + 1):
                                    self.writer.add_scalar(name + "_" + str(i), param.data.cpu().numpy().flatten()[i-1], batch_count)
                                    self.writer.add_scalar(name + "_grad_" + str(i), param.grad.cpu().numpy().flatten()[i-1], batch_count)
                            else:
                                self.writer.add_histogram(name, param.data.cpu().numpy().flatten(), batch_count)
                                self.writer.add_histogram(name + "grad", param.grad.cpu().numpy().flatten(), batch_count)
                                self.writer.flush()  ## remove this line, if results in slowdown, flush only at end

                # update batch count
                batch_count += 1
                epoch_loss += loss
                n += 1
            
            epoch = epoch+1
            # compute metrics
            if (epoch % eval_every == 0) and (self.evaluator is not None):
                metrics = {}
                metrics['train'] = self.evaluator.getMetrics(dataset,
                        eval_mode=True, report_threshold=True, threshold=0.5,
                        metrics_calculation=metrics_calculation, use_mask=use_mask
                )
                metrics['train']['loss'] = epoch_loss/(n + 1e-5)
                
                b = 0
                valid_loss = 0
                if validation_dataset is not None: 
                    for batch in validation_dataset:
                        oom = False
                        # up    date the model weights
                        batch_valid_data = processBatch(self.device, batch)
                        #batch, y, mask = batch_data['batch'], batch_data['y'], batch_data['mask']
                
                        # check for OOM errors
                        #loss = self.optimizer_step(batch_data, **optimizer_kwargs)
                        #try:
       

                        valid_loss += self.optimizer_step(batch_valid_data, epoch=epoch, valid=True, **optimizer_kwargs)
                        b+=1
                    metrics['validation'] = self.evaluator.getMetrics(validation_dataset,
                        eval_mode=True, threshold=0.5, metrics_calculation=metrics_calculation, use_mask=use_mask
                    )
                    
                    metrics['validation']['loss'] = valid_loss/(b + 1e-5)
                # report performance
                if not self.quiet:
                    reportMetrics(metrics,
                        label=epoch,
                        label_key='Epoch',
                        header=first_epoch
                    )
                self.updateHistory(metrics, epoch)
                first_epoch = False
                
                if best_state_metric:
                    state_metric = metrics[best_state_metric_dataset][best_state_metric]
                    if best_state_metric_goal == 'max' and state_metric > best_state_metric_threshold:
                        if state_metric > self.best_state_metric:
                            self.best_state_metric = state_metric
                            self.best_state = copy.deepcopy(self.model.state_dict())
                            self.best_epoch = epoch
                            self.metrics_history['best_epoch'] = epoch
                    elif best_state_metric_goal == 'min' and state_metric < best_state_metric_threshold:
                        if state_metric < self.best_state_metric:
                            self.best_state_metric = state_metric
                            self.best_state = copy.deepcopy(self.model.state_dict())
                            self.best_epoch = epoch
                            self.metrics_history['best_epoch'] = epoch
                
            # checkpoint
            if checkpoint_every and (epoch % checkpoint_every == 0):
                fname = self.saveState(epoch, "{}.{}.tar".format(self.model_name, epoch))
                logging.info("Writing checkpoint to file {} at epoch {}".format(fname, epoch))
            
        self.endTraining()
    
    def optimizer_step(self, batch_data, epoch, use_mask=True, use_weight=True, weight=None,
            valid=False):
        if not valid:
            self.optimizer.zero_grad()
        else:
            self.model.eval()
        output = self.model(batch_data['batch'])
        y = torch.cat((batch_data['y_pwm0'], batch_data['y_pwm1']), dim=0)
        y_mask = torch.cat((batch_data['pwm_mask0'], batch_data['pwm_mask1']), dim=0)
        
        out_mask = torch.cat((batch_data['dna_mask0'], batch_data['dna_mask1']), dim=0)

        
        ic_loss, weight = self.icLoss(output, y, out_mask, y_mask, background=[0.25,0.25,0.25,0.25])
        
        
        mse_loss = self.mseLoss(output, y, out_mask, y_mask, weight=self.rescaleWeight(weight),
                l1=True)
        
        seq = torch.cat((batch_data['y_hard0'], batch_data['y_hard1']), dim=0)
        seq_mask = torch.cat((batch_data['dna_mask0'], batch_data['dna_mask1']), dim=0)
        seq = seq[seq_mask]
        
        l = output.shape[0]//2
        #symm_loss = self.mseLoss(output[:l,:],output[l:,:].flip([0,1]), out_mask[:l],out_mask[:l],l1=True)
        
        #if epoch < 20:
        loss = self.mse_loss_weight*mse_loss + self.ic_loss_weight*ic_loss 
        #else:
        #    loss = self.mse_loss_weight*mse_loss + self.ic_loss_weight*ic_loss + symm_loss
        #+ self.seq_likelihood_loss(seq, output[out_mask], y[y_mask])
        
        if not valid:
            loss.backward()
            self.optimizer.step()
        else:
            self.model.train()
        return loss.item()
    
    def updateHistory(self, metrics, epoch):
        # update epoch
        self.metrics_history['epochs'].append(epoch)
        
        # update tags/metrics
        for tag in metrics:
            if tag not in self.metrics_history:
                self.metrics_history[tag] = {}
            
            for metric in metrics[tag]:
                # add metric to Tensorboard writer
                if(self.writer):
                    self.writer.add_scalar("{}/{}".format(tag, metric), metrics[tag][metric], epoch)
                
                # add metric to history
                if metric not in self.metrics_history[tag]:
                    self.metrics_history[tag][metric] = []
                self.metrics_history[tag][metric].append(metrics[tag][metric])
    
    def getHistory(self, tag, metric, epoch):
        if epoch == -1:
            ind = -1
        else:
            ind = self.metrics_history['epochs'].index(epoch)
        
        return self.metrics_history[tag][metric][ind]
    
    def saveState(self, epoch, suffix, metrics=True, optimizer=True, state=None):
        fname = ospj(self.checkpoint_path, suffix)
        
        # remove 'module' prefix from state_dict entries
        new_state_dict = OrderedDict()
        if state is None:
            state = self.model.state_dict()
        
        for k, v in state.items():
            name = k.replace("module.", "") # remove 'module.' prefix
            new_state_dict[name] = v
        
        data = {
            'model_state_dict': new_state_dict,
            'epoch': epoch
        }
        
        if metrics:
            data['history'] = self.metrics_history
            
        if optimizer:
            data['optimizer_state_dict'] = self.optimizer.state_dict()
        
        torch.save(data, fname)
        
        return fname
    
    def endTraining(self, message="Training Successfully Ended."):
        """Stuff we want to do at the end of training"""
        logging.info(message)
        
        # Save best state to file if we kept it
        if self.best_state is not None:
            fname = self.saveState(
                self.best_epoch, "{}.{}.tar".format(self.model_name, "best"),
                metrics=False,
                optimizer=False,
                state=self.best_state
            )
            logging.info("Writing best state to file {} (epoch: {})".format(fname, self.best_epoch))
            logging.info("Best tracked metric achieved: {:.3f}".format(self.best_state_metric))
        
        # Save metrics history to file
        logging.info("Saving metrics history to file.")
        logging.info("learned global temperature : {}".format(torch.sigmoid(self.model.global_temp).data.cpu().numpy()))
        MH = open(ospj(self.checkpoint_path, "{}_metrics.json".format(self.model_name)), "w")
        MH.write(json.dumps(self.metrics_history, indent=2))
        MH.close()

    def icLoss(self, output, y, out_mask, y_mask, background=[0.25,0.25,0.25,0.25]):
        background = torch.Tensor([background]*y[y_mask].shape[0]).to(self.device).log()
        eps = 1e-10
        term1 = self.kldiv(background, y[y_mask]).sum(dim=1)
        term2 = self.kldiv(background, torch.softmax(output[out_mask], dim=1)).sum(dim=1)
        
        ic_loss = torch.abs(term1 - term2).mean()
        return ic_loss, term1
    
    def mseLoss(self, output, y, out_mask, y_mask, weight=None, l1=False):
        softmax = torch.nn.functional.softmax
        if weight is None:
            if l1:
                return torch.nn.functional.l1_loss(softmax(output[out_mask], dim=1), y[y_mask])
            else:
                return torch.nn.functional.mse_loss(softmax(output[out_mask], dim=1), y[y_mask])
        else:
            if l1:
                loss = torch.nn.functional.l1_loss(softmax(output[out_mask], dim=1), y[y_mask], reduction='none').sum(dim=1)
            else:
                loss = torch.nn.functional.mse_loss(softmax(output[out_mask], dim=1), y[y_mask], reduction='none').sum(dim=1)
            return (weight * loss).mean()

    def ceLoss(self, output, y, out_mask, y_mask, weight=None):
        if weight is None:
            return self.criterion(output[out_mask], y[y_mask])
        else:
            loss = torch.nn.functional.cross_entropy(output[out_mask], y[y_mask], reduction='none')
            return (weight * loss).mean()

    def rescaleWeight(self, weight):
        ## rescale weight such that information contents in the middle gets most weight
        #weight = -(3/2)*(weight - (7/6))**2 + (49/24)
        #weight = (2 - (weight - 1)**2)
        weight = 1
        #weight.requires_grad = False
        return weight
    
    def seq_likelihood_loss(self, seq, out, true):
        #only take one strand
        eps = 1e-8
        l = len(seq)//2
        
        out_probs = (seq[:l,:]*torch.nn.functional.softmax(out[:l,:], dim =1)).sum(dim=1)
        true_probs = (seq[:l,:]*torch.nn.functional.softmax(true[:l,:], dim=1)).sum(dim=1)
        out_log = torch.log(out_probs + eps).sum()
        true_log = torch.log(true_probs + eps).sum()
        #print(out_log, true_log)

        diff_avg =  torch.abs(true_log - out_log)/l
        return diff_avg
