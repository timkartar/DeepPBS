import logging

def createFormattedStrings(fields, values, widths=None, pad=2, alignment='<', float_format='.3f'):
    """ automatically create a formatted string based on fields and values """
    if(widths is None):
        widths = [max(len(f)+pad, 6) for f in fields]
    header_format = ""
    values_format = ""
    for v, w in zip(values, widths):
        if(isinstance(v, int)):
            t = 'd'
        elif(isinstance(v, float)):
            t = float_format
        else:
            t = 's'
        header_format += '{:'+alignment+str(w)+'s}'
        values_format += '{:'+alignment+str(w)+t+'}'
    
    header_str = header_format.format(*fields)
    values_str = values_format.format(*values)
    
    return header_str, values_str, widths

def reportMetrics(metrics_dict, 
        label=None,
        label_key="label",
        label_width=None,
        sep_char=' | ',
        header=True,
        legend=True,
        header_sep=False,
        header_sep_char='-',
        logger='',
        **kwargs
    ):
    header_strs = []
    values_strs = []
    tags = []
    widths = []
    if label is not None:
        tags.append("")
        if label_width is None:
            l, v, w = createFormattedStrings([label_key], [label])
        else:
            l, v, w = createFormattedStrings([label_key], [label], [label_width])
        header_strs.append(l)
        values_strs.append(v)
        widths += w
    
    for tag, metrics in metrics_dict.items():
        keys = list(metrics.keys())
        keys.sort()
        values = [metrics[k] for k in keys]
        
        hs, vs, w = createFormattedStrings(keys, values, **kwargs)
        header_strs.append(hs)
        values_strs.append(vs)
        widths.append(sum(w))
        tags.append(tag)
    
    if header:
        if legend:
            # legend string
            ls = sep_char.join(["{:^{width}s}".format(t, width=w) for t, w in zip(tags, widths)])
            logging.getLogger(logger).info(ls)
            
        # header string
        hs = sep_char.join(header_strs)
        logging.getLogger(logger).info(hs)
        if header_sep:
            logging.getLogger(logger).info(header_sep_char*len(hs))
    
    fs = sep_char.join(values_strs)
    logging.getLogger(logger).info(fs)
