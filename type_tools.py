# Convert a string to the best type
def str_to_type(string):

    # Try integer
    try:    return int(string)
    except: pass

    # Try float
    try:    return float(string)
    except: pass

    # bool
    try:
        if string.lower().strip() == "true":  return True
        if string.lower().strip() == "false": return False
    except: pass
    
    # String
    return string
