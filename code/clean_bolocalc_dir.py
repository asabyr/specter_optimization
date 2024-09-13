import os
import glob
import sys

specter_dir=sys.argv[1]

first_layer=glob.glob(specter_dir+'/*.*', recursive=True)
second_layer=glob.glob(specter_dir+'/SPECTER/*.*', recursive=True)
third_layer=glob.glob(specter_dir+'/SPECTER/BF/*.*', recursive=True)
fourth_layer=glob.glob(specter_dir+'/SPECTER/BF/config/*.*', recursive=True)
fifth_layer=glob.glob(specter_dir+'/SPECTER/BF/config/Bands/Detectors/*.*', recursive=True)

for file_txt in first_layer:
    fname=file_txt.replace(specter_dir,'')
    
    if fname!='sensitivity.txt':
        os.remove(file_txt)

for file_txt in second_layer:
    fname=file_txt.replace(specter_dir+'/SPECTER/','')
    
    if fname!='sensitivity.txt':
        os.remove(file_txt)

for file_txt in third_layer:
    fname=file_txt.replace(specter_dir+'/SPECTER/BF/','')
    
    if fname!='sensitivity.txt' and fname!='optical_power.txt':
        os.remove(file_txt)

for file_txt in fourth_layer:
    fname=file_txt.replace(specter_dir+'/SPECTER/BF/config/','')
    
    if fname!='camera.txt' and fname!='channels.txt' and fname!='optics.txt':        
        os.remove(file_txt)

for file_txt in fifth_layer:
    os.remove(file_txt)
