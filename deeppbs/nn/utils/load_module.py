import importlib.util
import os.path as op

def loadModule(model_name, module_path, model_args, model_kwargs):
    # specify the module that needs to be 
    # imported relative to the path of the 
    # module
    name = op.basename(module_path).strip('.py')
    spec = importlib.util.spec_from_file_location(name, module_path)
    
    # creates a new module based on spec
    module = importlib.util.module_from_spec(spec)
    
    # executes the module in its own namespace
    # when a module is imported or reloaded.
    spec.loader.exec_module(module)
    
    # initialize the model
    model = getattr(module, model_name)(*model_args, **model_kwargs)
    
    return module, model
