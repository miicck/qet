open_logfiles = {}
def log(message, filename="qet.out"):

    # Open the log file in line-buffered mode
    if not filename in open_logfiles:
        open_logfiles[filename] = open(filename, "w", 1)

    # Write the message
    open_logfiles[filename].write(str(message) + "\n")
