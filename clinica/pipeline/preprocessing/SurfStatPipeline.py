#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 15:20:40 2016

@author: Junhao WEN
"""

def SurfStatGruopAnalysis(ContrastLinearModel,Format, CSVFilename, PATH_TO_RECON_ALL_OUTPUTS, a_required_path, result_dir):
    """
        This is to use SurfStat to do the Group analysis for the reconAll outputs, after the reconAll pipeline, you should just define the paths to 
        SurfStatGruopAnalysis, and create the CSV file, and run the pipeline, at last, you will get the results images.
        
        Inputnode
        ---------
        SurfStat
        Inputs: ContrastLinearModel: string, the linear model that fit into the GLM, for example '1+Lable'.
                Format: string, the format that you want to use for your CSV file column variables, it depends on the CSV file. 
                CSVFilename: string, the path to your csv file.
                PATH_TO_RECON_ALL_OUTPUTS:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh.
                a_required_path: this is the path to find the matlabscript(NipypeSurfStat.m) that you are going to run.
                result_dir: the directory to contain the result images.
                
          For more infomation about SurfStat, please check:
          http://www.math.mcgill.ca/keith/surfstat/
         
      Outputnode:
      
      :param: ContrastLinearModel: string, the linear model that fit into the GLM, for example '1+Lable'.
      :param: Format: string, the format that you want to use for your CSV file column variables, it depends on the CSV file. 
      :param: CSVFilename: string, the path to your csv file.
      :param: PATH_TO_RECON_ALL_OUTPUTS:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh.
      :param: a_required_path:  this is the path to find the matlabscript(NipypeSurfStat.m) that you are going to run.
      :param: result_dir: the directory to contain the result images.
      
      return images in result_dir
    
    """
    from nipype.interfaces.utility import Function
    import nipype.pipeline.engine as pe 
     
    def runmatlab(ContrastLinearModel, Format, CSVFilename, PATH_TO_RECON_ALL_OUTPUTS, a_required_path, result_dir ):
        from nipype.interfaces.matlab import MatlabCommand, get_matlab_command
        MatlabCommand.set_default_matlab_cmd(get_matlab_command())#this is to set the matlab_path in your bashrc file, if the user r not in bash or linux, but anyway, this is just to find ur matlab
        matlab = MatlabCommand()
        matlab.inputs.paths = [a_required_path]# this is the path to add into matlab, addpath
        matlab.inputs.script = """
        NipypeSurfStat('%s', '%s', '%s', '%s', '%s', '%s')"""%(ContrastLinearModel, Format, 
                                     CSVFilename, PATH_TO_RECON_ALL_OUTPUTS, a_required_path, result_dir )  # here, we should define the inputs for the matlab function that you want to use
        matlab.inputs.mfile = True # this will create a file: pyscript.m , the pyscript.m is the default name
        out = matlab.run()
        return out
    
    SurfStat = pe.Node(name='SurfStat',
                   interface=Function(input_names=['ContrastLinearModel', 'Format', 'CSVFilename', 
                                         'PATH_TO_RECON_ALL_OUTPUTS', 'a_required_path', 'result_dir'],
                                      output_names=['out_file'],
                                      function=runmatlab))
    SurfStat.inputs.ContrastLinearModel = ContrastLinearModel
    SurfStat.inputs.Format = Format
    SurfStat.inputs.CSVFilename = CSVFilename
    SurfStat.inputs.PATH_TO_RECON_ALL_OUTPUTS = PATH_TO_RECON_ALL_OUTPUTS
    SurfStat.inputs.a_required_path = a_required_path
    SurfStat.inputs.result_dir = result_dir
    
    return SurfStat
