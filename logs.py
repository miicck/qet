import datetime

open_logfiles = {}
def log(message, filename="qet.out"):

    # Open the log file in line-buffered mode
    if not filename in open_logfiles:
        new_open = open(filename, "w", 1)
        open_logfiles[filename] = new_open
        new_open.write("\nLogfile (re)opened: "+str(datetime.datetime.now())+"\n")

    # Write the message
    open_logfiles[filename].write(str(message) + "\n")
