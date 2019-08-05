import time
from IPython.display import clear_output

class Timer():
    
    def start(self):
        self.t0 = time.time()
        #print(self.t0)
        
    def stop(self, message):
        self.t1 = time.time()
        print(message, self.t1 - self.t0, 'sec')
        
def progress_bar(progress):
    bar_length = 20
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1
        
    block = int(round(bar_length * progress))
    clear_output(wait = True)
    text = "Progress: [{0}] {1:.1f}%".format( "#" * block + "-" * (bar_length - block), progress * 100)
    print(text)
        
       