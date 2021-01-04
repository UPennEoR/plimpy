#!/lustre/aoc/projects/hera/jaguirre/anaconda3/envs/hera/bin/python

import sys, getopt, os

def main(argv):
   JD = ''
   try:
      opts, args = getopt.getopt(argv,"hj:",["jd="])
   except getopt.GetoptError:
      print('test.py -j <julian_date>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('runLSTSliceNight.py -j <julian_date>')
         sys.exit()
      elif opt in ("-j", "--jd"):
         JD = arg
   print('Julian date is %s' % JD)

   if JD != '':

       import LSTSliceNight
       LSTSliceNight.sliceJD(JD)

if __name__ == "__main__":
   main(sys.argv[1:])
