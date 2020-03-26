import datetime

# Log a message in a logfile with the given name
# these logfiles are at a persistant location and
# are created in whatever the working directory 
# was when the file was first logged to.
open_logfiles = {}
def log(message, filename="qet.out"):
    global open_logfiles

    # Open the log file in line-buffered mode
    if not filename in open_logfiles:
        new_open = open(filename, "w", 1)
        open_logfiles[filename] = new_open
        new_open.write("\nLogfile opened: "+str(datetime.datetime.now())+"\n")

    # Write the message
    open_logfiles[filename].write(str(message) + "\n")
