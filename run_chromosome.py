#!/usr/bin/env python
import socket
import os
import shutil
FNULL = open(os.devnull, 'w')
SPARX = 'sparx-'+socket.gethostname()[0:2]
PRESPARX = 'presparx-'+socket.gethostname()[0:2]
WORKING_DIR='/tiara/home/ithsieh/N1333_GA/working_dir'


from subprocess import call

def run( par, par_fixed, chromosome_id, target):
    pid_str = str(chromosome_id)
    DIR_ID='chromosome_'+pid_str
    
    call(['mkdir','-p','working_dir/'+DIR_ID])
    with open('preprocessor/model.py', 'rt') as fin:
        with open("working_dir/"+DIR_ID+"/model.py", "wt") as fout:
            for line in fin:
                new_line = line
                
                for key in par.keys():
                    new_line = new_line.replace(str(key), str(par[key]) )
                
                for key in par_fixed.keys():
                    new_line = new_line.replace(str(key), str(par_fixed[key]) )
                
                fout.write(new_line)

    # enter working directory
    os.chdir(WORKING_DIR+'/'+DIR_ID)

    # generate model
    
    shutil.copy('../../preprocessor/grid.py','grid.py')
    call([PRESPARX,'-o','model'], stdout=FNULL)
    
    
    result = []
    # loop for wavelength
    for i in range(len(target.wavelen)):
        image_file='cont_map_%g'%(target.freq[i])
        call([SPARX,'run','task_contobs',
              'source=model',
              'out='+image_file,
              'wavelen='+str(target.wavelen[i])+'m',
              'dist='+target.distance,
              "cell=['"+target.cell+"','"+target.cell+"']",
              'npix=['+str(target.npix)+','+str(target.npix)+']',
              "unit=JY/PIXEL",
              'subres='+target.subres]
             , stdout=FNULL 
            )
        # loop fot beamsize
        result.append( [] )
        
        for j in range(len(target.luminance[i])):
            major_size = target.luminance[i][j][0]
            minor_size = target.luminance[i][j][1]
            convolved_file='conv_map%g_%g_%g' % ( target.freq[i], major_size, minor_size)
            call(['convol',
                  'map='+image_file,
                  'out='+convolved_file,
                  'fwhm=%g,%g'%(major_size,minor_size)]
                 , stdout=FNULL 
                )
            call(['imspect',
                  'in='+convolved_file,
                  "region=rel,box(0,0,0,0)",
                  'log=tmp.log']
                 , stdout=FNULL 
                )
            with open('tmp.log', 'r') as f:
                for line in f:
                    pass
                last = line
            simulated_luminance = float(last.rstrip('\n').split()[2])
            result[i].append( simulated_luminance )
            
            #os.remove('tmp.log')
            call(['rm','-rf','tmp.log'])
            
            #shutil.rmtree(convolved_file)
            call(['rm','-rf',convolved_file])
        #shutil.rmtree(image_file)
        #os.remove(image_file+'.fits')
        call(['rm','-rf',image_file,image_file+'.fits'])
    #os.remove('model')
    call(['rm','-rf','model'])

    

    # leave the directory
    os.chdir('../..')
    
    return result
