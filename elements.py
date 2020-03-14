import os
from qet.types import str_to_type 

# Read data on elements from elements.csv
def get_elements():

    data = {}
    headings = None

    path = os.path.dirname(__file__)
    with open(path+"/elements.csv") as f:
        for line in f:

            if headings is None:
                headings = line.split(",")
                continue

            elm = [str_to_type(w) for w in line.split(",")]
            if len(elm) > len(headings):
                elm = elm[0:len(headings)]

            # Data is of the format
            # data[symbol][quantity]
            data[elm[2]] = {}
            for i in range(0, len(elm)):
                data[elm[2]][headings[i]] = elm[i]

    return data

elements = get_elements()
