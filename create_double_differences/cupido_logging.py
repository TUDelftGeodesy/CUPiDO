## cupido_logging.py
#
#   This function handles the logging of the cupido function.
# 
#   (c) Sami Samiei Esfahany, Adriaan van Natijne, Hans van der Marel and Freek van Leijen
#       Delft University of Technology, 2016.
#
#   Version:    1.0 
#   Created:    19 September 2016
#   Modified:   
#



import time

class logging(object):
    def __init__(self, base_name, extention, active):
        self.active = active
        if(active):
            self.time_stamp = time.localtime()
            self.file_name = base_name + '_' + time.strftime("%a, %d %b %Y %H%M%S +0000",self.time_stamp).replace(" ", "_").replace(",", "").replace("+", "") + '.' + extention
            self.file_handle = open(self.file_name, 'w')
            self.write('CUPiDO logfile created at ' + time.strftime("%a, %d %b %Y %H:%M:%S +0000",self.time_stamp))

    def write(self, message):
        if(self.active):
            self.file_handle.write(message + '\n')
        print(message)

