import os, glob

os.system("rm -fr  ../PlotterSteer/*.steer")
os.system("cp /nfs/dust/cms/user/"+os.environ["USER"]+"/WorkingArea/SFramePlotter/Analysis/* ../PlotterSteer/")
for f in glob.glob("../PlotterSteer/*.steer"):
    os.system("mv "+f+" "+f.replace(".steer","_store.steer"))
os.system("ln -s /nfs/dust/cms/user/"+os.environ["USER"]+"/WorkingArea/SFramePlotter/Analysis/* ../PlotterSteer/")
