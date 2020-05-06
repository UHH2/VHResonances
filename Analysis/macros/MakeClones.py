import os

os.system("rm -fr  ../Plotter/*.steer")
os.system("ln -s /nfs/dust/cms/user/"+os.environ["USER"]+"/WorkingArea/SFramePlotter/Analysis/* ../Plotter/")
