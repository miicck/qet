# Convert a string to the best type
def str_to_type(string):

    # integer or float
    try:
        fval = float(string)
        ival = int(string)
        # If the floating point value differs
        # this was actually a float that has 
        # been truncated
        if fval != float(ival):
            # A valid float
            return fval
        else:
            # A valid integer
            return ival

    except: pass

    # bool
    try:
        if string.lower().strip() == "true":
            return True

        if string.lower().strip() == "false":
            return False
    except: pass
    
    # String
    return string
