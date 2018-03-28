def debug_breakpt(option='ipdb'):
    if option == 'ipdb':
        import ipdb
        ipdb.set_trace()
    else:
        raise NotImplementedError