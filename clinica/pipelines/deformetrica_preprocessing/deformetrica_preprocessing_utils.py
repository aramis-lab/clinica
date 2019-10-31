# -*- coding: utf-8 -*-

__author__ = "Alexis Guyot"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Alexis Guyot"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Alexis Guyot"
__email__ = "alexis.guyot@icm-institute.org"
__status__ = "Development"



def select_object(
        participant_id,
        session_id,
        caps_layout,
        object_name,
        object_description):
    """Select subject-specific object

    Returns user-specified object corresponding to participant ID + session ID

    Args:
       participant_id (String):
       session_id (String)
       caps_layout (CAPS layout object):
       object_name (String): Freesurfer identifier for the object
           (e.g. 'aseg', 'brainmask')
       object_description (String): user-specified object description
           (e.g. 'segmentation', 'brain mask')

    Returns:
       object_path (String): full path to the object identified by
           participant ID + session ID
    """
    import warnings
    import os

    # make participant_id compliant with caps_layout: crop 'sub' and
    # 'ses'
    if participant_id.startswith('sub-'):
        participant_id = participant_id[4:]
    if session_id.startswith('ses-'):
        session_id = session_id[4:]

    # search for object with the caps_layout
    selected_objects_unfiltered = caps_layout.get(
        return_type='file',
        freesurfer_file=object_name,
        extensions='\\.mgz',
        subject=participant_id,
        session=session_id,
        regex_search=True)

    # this will return several objects eg. aparc+[object_name].mgz,
    # [object_name].auto.mgz, [object_name].mgz ... We only want
    # '[object_name].mgz', hence a manual filter
    selected_objects = []
    for obj in selected_objects_unfiltered:
        object_filename = os.path.basename(obj)
        if object_filename == '{0}.mgz'.format(object_name):
            selected_objects.append(obj)

    # check for errors/abnormalities
    if len(selected_objects) == 0:
        raise RuntimeError(
            'No {0}s were found for participant {1} and session {2}'.format(
                object_description, participant_id, session_id))
    if len(selected_objects) > 1:
        warnings.warn(
            'Several {0}s were found for participant {1} and session {2}'.format(
                object_description, participant_id, session_id),
            RuntimeWarning)

    return selected_objects[0]



def select_caps_brains(subjects, sessions, caps_layout):
    """Returns all brain mask volumes for a list of subjects/sessions

    Args:
       subjects (list of String)
       sessions (list of String)
       caps_layout (CAPS layout object)

    Returns:
       brain_list (list of String): list of all the .mgz brain masks
    """
    if len(subjects) != len(sessions):
        raise RuntimeError(
            "Subjects list and sessions list must have the same length.")
    brain_list = [
        select_brain(
            subjects[i],
            sessions[i],
            caps_layout) for i in range(len(subjects))]

    return brain_list


def select_brain(participant_id, session_id, caps_layout):
    """Select subject-specific brainmask.mgz

    Returns brainmask.mgz corresponding to a participant ID + session ID

    Args:
       participant_id (String):
       session_id (String)
       caps_layout (CAPS layout object):

    Returns:
       brain (String): brainmask.mgz segmentations full path
    """
    brain = select_object(
        participant_id,
        session_id,
        caps_layout,
        'brainmask',
        'brain mask')

    return brain



def select_caps_segmentations(subjects, sessions, caps_layout):
    """Select segmentation for list of subjects

    Returns all segmentation volumes for a list of subjects/sessions

    Args:
       subjects (list of String)
       sessions (list of String)
       caps_layout (CAPS layout object)

    Returns:
       segmentation_list (list of String): list of all the aseg.mgz
           segmentations
    """
    if len(subjects) != len(sessions):
        raise RuntimeError(
            "Subjects list and sessions list must have the same length.")
    segmentation_list = [
        select_segmentation(
            subjects[i],
            sessions[i],
            caps_layout) for i in range(len(subjects))]

    return segmentation_list



def select_segmentation(participant_id, session_id, caps_layout):
    """ Select subject-specific aseg.mgz

    Returns aseg.mgz segmentation corresponding to a
    participant ID + session ID

    Args:
       participant_id (String):
       session_id (String)
       caps_layout (CAPS layout object):

    Returns:
       segmentation (String): aseg.mgz segmentations full path
    """
    segmentation = select_object(
        participant_id,
        session_id,
        caps_layout,
        'aseg',
        'segmentation')

    return segmentation



def get_file(filename):
    """This only checks whether the file exists

    Args:
        filename (String): path to be checked

    Returns:
        filename (String): filename (if file exists)
    """
    import os.path

    if not os.path.isfile(filename):
        raise RuntimeError("{0} is not a valid file".format(filename))

    return filename



def select_colin27_image(colin27_image):
    """Loads Colin 27 volume

    Args:
        colin27_image (String): path to Colin27 MR volume

    Returns:
        colin27_image (String): path to Colin27 MR volume (if exists)
    """

    return get_file(colin27_image)



def check_cuda_flag(in_cuda_flag):
    """ Check the CUDA flag provided by the user

    Args:
        in_cuda_flag (String): CUDA flag provided by the user

    Returns:
        out_cuda_flag (Boolean): yes of in_cuda_flag in
            ['Y', 'y', 'YES', 'yes', 'Yes', '1', 'True', 'TRUE', 'true']
    """
    out_cuda_flag = False
    yes_list = ['Y', 'y', 'YES', 'yes', 'Yes', '1', 'True', 'TRUE', 'true']
    if in_cuda_flag in yes_list:
        out_cuda_flag = True

    return out_cuda_flag



def read_surface_points(in_surface):
    """
    Reads points from a vtk surface

    Args:
        in_surface (vtk surface)

    Returns:
        out_points (numpy array): [n,3] array
    """
    import numpy as np

    point_number = in_surface.GetNumberOfPoints()
    out_points = np.zeros((point_number, 3))
    for point_index in range(point_number):
        point_x, point_y, point_z = in_surface.GetPoint(point_index)
        out_points[point_index, 0] = point_x
        out_points[point_index, 1] = point_y
        out_points[point_index, 2] = point_z

    return out_points



def get_target_coordinates(
        in_points,
        in_source_filename,
        in_target_filename,
        in_affine_filename,
        in_mm=True):
    """
    Function to get the coordinates of a set of points displaced with
    an affine transformation found using FSL flirt. This calls
    img2imgcoord, an FSL utility.

    Args:
        in_points (Numpy array): [n1, 3] array of points
        in_source_filename (String): path to the source volume
        in_target_filename (String): path to the target volume
        in_affine_filename (String): path to a .mat file containing a
            4x4 affine matrix
        in_mm (Boolean): flag indicating if going from [voxel to voxel]
            or [mm to mm] positions

    Returns:
        out_points (Numpy array): [n1, 3] array of points
    """
    import os
    import numpy as np
    # save input points to file
    in_points_filename = "./temp_inpoints_img2imgcoord.txt"
    np.savetxt(in_points_filename, in_points, fmt="%f")
    # apply img2imgcoord
    out_points_filename = "./temp_outpoints_img2imgcoord.txt"
    mm_flag = ""
    if in_mm:
        mm_flag = " -mm"
    img2imgcoord_command = "img2imgcoord -src {0} -dest {1} -xfm {2} {3}{4} > {5}".format(
        in_source_filename,
        in_target_filename,
        in_affine_filename,
        in_points_filename,
        mm_flag,
        out_points_filename)
    os.system(img2imgcoord_command)
    # read img2imgcoord output, skip first line
    with open(out_points_filename, 'r') as out_points_file:
        lines = (line for line in out_points_file if not line.startswith('C'))
        out_points = np.loadtxt(lines)

    return out_points



def check_structure_in_dictionary(structure, dictionary):
    """Check if a given structure is part of a dictionary.

    Raise an error if not.

    Args:
        structure (String): name of the structure to be made a template
            of
        dictionary (Dictionary): name of the dictionary of possible
            structures

    Returns:
        N/A
    """
    structure_list = dictionary.keys()
    structure_list_string = ", ".join(structure_list)
    if structure not in structure_list:
        raise RuntimeError(
            '\'{0}\' is not a supported structure.\n Valid structures are: {1}'.format(
                structure, structure_list_string))



def get_structure_id_dictionary():
    """Get structure Freesurfer ID
    
    Returns the dictionary that maps each hippocampal structure to a 
    freesurfer-intelligible ID.
    We may want to move this out of code and into file at a later stage

    Args:
        N/A

    Returns:
        structure_id_dictionary (dictionary): dictionary mapping each
            hippocampal structure to its corresponding freesurfer ID
    """
    structure_id_dictionary = {}
    # ventricle
    structure_id_dictionary["left-ventricle"] = 4
    structure_id_dictionary["right-ventricle"] = 43
    # thalamus
    structure_id_dictionary["left-thalamus"] = 10
    structure_id_dictionary["right-thalamus"] = 49
    # caudate
    structure_id_dictionary["left-caudate"] = 11
    structure_id_dictionary["right-caudate"] = 50
    # putamen
    structure_id_dictionary["left-putamen"] = 12
    structure_id_dictionary["right-putamen"] = 51
    # pallidum
    structure_id_dictionary["left-pallidum"] = 13
    structure_id_dictionary["right-pallidum"] = 52
    # hippocampus
    structure_id_dictionary["left-hippocampus"] = 17
    structure_id_dictionary["right-hippocampus"] = 53
    # amygdala
    structure_id_dictionary["left-amygdala"] = 18
    structure_id_dictionary["right-amygdala"] = 54

    return structure_id_dictionary



def get_colin27_structure_path_dictionary(in_colin27_resources_folder):
    """Get map Colin 27 hippocampal structure -> path

    Returns the dictionary that maps each hippocampal structure to the
    path of the corresponding structure for Colin 27.
    We may want to move this out of code and into file at a later stage

    Args:
        in_colin27_resources_folder (String): path to folder with
            Colin 27 resources

    Returns:
        colin27_structure_path_dictionary (dictionary): dictionary
            mapping each hippocampal structures to the path of their
            corresponding Colin 27 structure
    """
    colin27_structure_path_dictionary = {}
    side_list = [
        'left',
        'right']
    structure_list = [
        'ventricle',
        'thalamus',
        'caudate',
        'putamen',
        'pallidum',
        'hippocampus',
        'amygdala']
    for structure in structure_list:
        for side in side_list:
            side_structure = '{0}-{1}'.format(side, structure)
            colin27_structure_path_dictionary[side_structure] = "{0}/colin27_{1}.vtk".format(
                in_colin27_resources_folder,
                side_structure)

    return colin27_structure_path_dictionary



def get_valid_structure_list(in_colin27_resources_folder):
    """Returns the list of structures that are currently supported.

    Takes the intersection of:
        - structures for which a freesurfer ID has been defined
        - structures for which a path to corresponding Colin 27
            structures has been defined

    Args:
        in_colin27_resources_folder (String): path to folder with
            Colin 27 resources

    Returns:
        valid_structure_list (list of string): list of all valid
            structures
    """
    import deformetrica_preprocessing_utils as utils

    # get structures for which a freesurfer ID has been defined
    freesurferid_structure_dictionary = utils.get_structure_id_dictionary()
    freesurferid_structure_list = freesurferid_structure_dictionary.keys()
    freesurferid_structure_set = set(freesurferid_structure_list)

    # get structures for which a path to corresponding Colin 27
    # structures has been defined
    colin27_structure_path_dictionary = utils.get_colin27_structure_path_dictionary(
        in_colin27_resources_folder)
    colin27_structure_path_list = colin27_structure_path_dictionary.keys()
    colin27_structure_path_set = set(colin27_structure_path_list)

    # intersect the two sets
    valid_structure_set = colin27_structure_path_set.intersection(
        freesurferid_structure_set)
    valid_structure_list = list(valid_structure_set)

    return valid_structure_list



def get_structure_id(in_structure):
    """Returns a freesurfer-intelligible ID for a given structure.

    All potential IDs are stored in a dictionary

    Args:
        in_structure (String): Structure name

    Returns:
        out_structure_id (integer): corresponding freesurfer ID (>=0)
    """
    import deformetrica_preprocessing_utils as utils

    # get dictionaries linking structures to Freesurfer ID
    structure_id_dictionary = utils.get_structure_id_dictionary()
    # check if structure is valid
    utils.check_structure_in_dictionary(in_structure, structure_id_dictionary)
    # get structure ID
    out_structure_id = structure_id_dictionary[in_structure]

    return out_structure_id



def get_colin27_structure_template(in_structure, in_colin27_resources_folder):
    """ Get path to Colin 27 structure template
    
    Returns for each structure ID a path to a corresponding Colin27
    template.
    This essentialy returns a path to Clinica resource folder.
    To do: we may want to do some processing on the .vtk contained in
           Clinica resource folder (e.g. extra smoothing) and store it
           inside the node mapflow in the future.
    All potential paths are stored in a dictionary

    Args:
        in_structure (String): Structure name
        in_colin27_resources_folder (String): path to folder with
            Colin27 resources

    Returns:
        out_colin27_structure_template (String): path to corresponding
            Colin27 template (for now: path to corresponding Colin27
            structure in resource folder)
    """
    import deformetrica_preprocessing_utils as utils

    # define dictionaries linking structures to Colin 27 counterpart
    # We may want to move this out of code and into file in the future
    colin27_structure_path_dictionary = utils.get_colin27_structure_path_dictionary(
        in_colin27_resources_folder)
    # check if structure is valid
    utils.check_structure_in_dictionary(
        in_structure,
        colin27_structure_path_dictionary)
    # get structure ID
    out_colin27_structure_path = utils.get_file(
        colin27_structure_path_dictionary[in_structure])
    # for now: no extra data manipulation
    out_colin27_structure_template = out_colin27_structure_path

    return out_colin27_structure_template



def process_structure(
        in_structure,
        in_structure_file,
        in_colin27_resources_folder):
    """Returns a list (with potentially a single element) of structures

    Args:
        in_structure (String): name of the structure to be made a
            template of
        in_structure_file (String): path to file containing a list of
            names of the structures to be made a template of
        in_colin27_resources_folder (String): path to folder with
            Colin 27 resources

    Returns:
        out_structure_list (list of string): list of the structures to
            be processed. Will only contain one element if the structure
            was passed with flag '-st'/'--structure'. Will have 1+
            (parsed) elements if a file was provided with flag
            '-stf' / '--structure-file'
        out_structure_file (String): File containing the list of
            structures. Structures are listed line by line.
    """
    import os
    import deformetrica_preprocessing_utils as utils

    # retrieve list of structures from the arguments
    out_structure_list = []
    if in_structure is None and in_structure_file is None:
        valid_structure_list = utils.get_valid_structure_list(
            in_colin27_resources_folder)
        valid_structure_list_string = ", ".join(valid_structure_list)
        raise RuntimeError(
            'Please provide a structure in: {1}'.format(
                valid_structure_list_string))
    elif in_structure_file is not None:
        # This takes precedence over in_structure, ie if both flags are
        # provided, only 'stf'/'--structure-file' will be taken into
        # account
        # parse structure file
        with open(in_structure_file, 'r') as structure_f:
            # populate output
            out_structure_list = structure_f.read().split()
    elif in_structure is not None:
        # populate output
        out_structure_list = [in_structure]

    # Write list of structure to file
    # Note: this 'kind of' returns the same file as in_structure_file
    #       (if provided), except that the formatting of
    #       in_structure_file can vary, while this of out_structure_file
    #       will always be the same (each structure written on a
    #       separate line)
    current_path = os.path.abspath('./')
    out_structure_file = '{0}/out_structure_file.txt'.format(current_path)
    with open(out_structure_file, 'w') as structure_f:
        for out_structure in out_structure_list:
            structure_f.write("{0}\n".format(out_structure))

    return [out_structure_list, out_structure_file]



def link_objects1_objects2(in_object1_list, in_object2_list, output_map=False):
    """ Links two lists of objects together.

    Links a list of objects (type 1) to another list of objects
    (type 2). Returns all possible combinations (Cartesian product)
    between the two lists.
    e.g. if in_object1_list = (A1, A2, A3) and
    in_object2_list = (B1, B2), then the function returns:
    out_object1_list    out_object2_list
    ----------------    ----------------
    A1                  B1
    A1                  B2
    A2                  B1
    A2                  B2
    A3                  B1
    A3                  B2
    The goal of this function is to allow the use of a MapNode over two
    separate iterables [object1] and [object2] by using a MapNode over
    the iterable [(object1, object2)]
    The function can also return a list of indices mapping each element
    of the input objects to a list of corresponding (replicated)
    elements in the output objects. This can be used later to 'unlink'
    the two objects, e.g.  for a given subject, get all corresponding
    structures / associated meshes, etc. after some processing carried
    out on a MapNode over the iterable [(subject, structure)]

    Args:
        in_object1_list (list): list of objects 1
        in_object2_list (list): list of objects 2
        get_unlink (Boolean): Flag to return the a list of indices
            mapping each element of in_object1_list to a list of
            elements of out_object1_list

    Returns:
        out_object1_list (list): reordered list of objects 1
        out_object2_list (list): reordered list of objects 2
        out_map (list): list that maps each element of in_object1_list
            to a list of elements of out_object1_list.
    """
    import itertools
    import numpy as np

    # generate all possible combinations of (object1, object2)
    object1_object2_list = list(
        itertools.product(in_object1_list, in_object2_list))
    out_object1_list = [
        object1_object2[0] for object1_object2 in object1_object2_list]
    out_object2_list = [
        object1_object2[1] for object1_object2 in object1_object2_list]

    # retrieve indices mapping objects from in_object1_list to
    # corresponding objects in out_object1_list
    out_map = None
    if output_map:
        out_object1_np = np.array(out_object1_list)
        out_map = [
            np.where(out_object1_np == in_object1)[0].tolist() for in_object1 in in_object1_list]

    return out_object1_list, out_object2_list, out_map




def link_objects_structures(
        in_brain_list,
        in_segmentation_list,
        in_affine_list,
        in_structure_list):
    """ Links object list to structure ID list

    Links a list of objects to a list of structures IDs. The list of
    object can be:
    - a list of path to brain images
    - a list of path to brain segmentations
    - a list of affine transformations

    Args:
        in_segmentation_list (list of String): list of paths to
            segmentations
        in_affine_list (list of String): list of paths to affine reg
            matrix
        in_brain_list (list of String): list of paths to the subject
            brain images
        in_structure_list (list of String): list of structures

    Returns:
        out_brain_list (list of String): reordered list of paths to the
            subject brain images
        out_segmentation_list (list of String): reordered list of paths
            to segmentation
        out_affine_list (list of String): reordered list of paths to
            affine reg matrix
        out_structure_list (list of String): reordered list of
            structures
        out_map (list of list of integers): list mapping each input
            subject to corresponding output subjects
    """
    import numpy as np
    import deformetrica_preprocessing_utils as utils

    # link [brain image] to [structure]
    [out_brain_list, out_structure_list, out_map] = utils.link_objects1_objects2(
        in_brain_list,
        in_structure_list,
        output_map=True)
    # link [segmentation] to [structure]
    [out_segmentation_list, dummy, dummy] = utils.link_objects1_objects2(
        in_segmentation_list,
        in_structure_list)
    # link [affine matrix] to [structure]
    [out_affine_list, dummy, dummy] = utils.link_objects1_objects2(
        in_affine_list,
        in_structure_list)

    return [out_brain_list, out_segmentation_list, out_affine_list, out_structure_list, out_map]



def apply_affine(
        in_segmentation,
        in_affine,
        in_brain,
        in_colin27,
        in_structure):
    """Affine reorient segmentation

    Apply an affine matrix (usually found with an affine registration) to
    reorient a segmentation.

    Args:
        in_segmentation (String): path to .vtk segmentation
        in_affine (String): path to text file containing 4x4 matrix
        in_brain (String): path to brain volume
        in_colin27 (String): path to colin27 volume
        in_structure (String): name of the structure that corresponds
            to the segmentation. This is not used to carry out the
            registration.  Instead it allows us to insert the structure
            name in the name of the file which is returned. This
            structure name will later be picked in a datasink regular
            expression.

    Returns:
        out_affine_reoriented (String): path to .vtk re-oriented segmentation
    """
    import os
    import vtk
    import numpy as np
    import deformetrica_preprocessing_utils as utils

    # Read input segmentation
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(in_segmentation)
    reader.Update()
    segmentation = reader.GetOutput()

    # Read .vtk points
    initial_points_RAS = utils.read_surface_points(segmentation)
    point_number = initial_points_RAS.shape[0]

    # go from T1 subject coordinate system to Colin27 coordinate system
    # using affine matrix
    transformed_points_RAS = utils.get_target_coordinates(
        initial_points_RAS,
        in_brain,
        in_colin27,
        in_affine,
        in_mm=True)

    # update .vtk surface with new points
    current_path = os.path.abspath('./')
    out_affine_reoriented = "{0}/affine_registered_{1}.vtk".format(
        current_path,
        in_structure)
    vtk_points = vtk.vtkPoints()
    for point_index in range(point_number):
        vtk_points.SetNumberOfPoints(point_number)
        vtk_points.SetPoint(
            point_index,
            transformed_points_RAS[point_index, 0],
            transformed_points_RAS[point_index, 1],
            transformed_points_RAS[point_index, 2])
    segmentation.SetPoints(vtk_points)
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(out_affine_reoriented)
    writer.SetInputData(segmentation)
    writer.Update()

    return out_affine_reoriented



def decimate_mesh(in_mesh, in_decimation):
    """ Decimate mesh

    Args:
        in_mesh (String): path to .vtk mesh
        in_decimation (float): decimation percentage (between 0 and 1)

    Returns:
        out_decimated_mesh (String): path to .vtk decimated mesh
    """
    import vtk
    import os

    # Read input
    # segmentation
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(in_mesh)
    reader.Update()

    # create decimate filter
    decimate=vtk.vtkDecimatePro()
    decimate.SetInputConnection(reader.GetOutputPort())
    decimate.SetTargetReduction(in_decimation)
    decimate.Update()

    # Write decimated surface to file
    current_path = os.path.abspath('./')
    out_decimated_mesh = '{0}/out_decimated_mesh.vtk'.format(current_path)
    writer = vtk.vtkPolyDataWriter()
    writer.SetInputConnection(decimate.GetOutputPort())
    writer.SetFileName(out_decimated_mesh)
    writer.Update()

    return out_decimated_mesh



def get_persubject_meshes(in_mesh_list, in_map):
    """ Returns per-subject list of decimated meshes.

    For each mesh in a list, will work out which subject it is
    associated to, so meshes extracted from the same subject can be
    grouped prior to rigid registration to Colin 27.

    Args:
        in_mesh_list (list of .vtk meshes): Decimated meshes that will
            be input into the rigid registration to Colin 27 algorithm
        in_map (list of list of integers): list mapping each input
            subjects to corresponding subjects in the list
            (input_subject, structure) - equivalent to (decimated mesh)

    Returns:
        out_persubject_meshlist (list of list of .vtk meshes):
            per-subject list of decimated meshes
    """
    out_persubject_meshlist = [
        [in_mesh_list[
            subject_index] for subject_index in subject_indices] for subject_indices in in_map]

    return out_persubject_meshlist



def get_mesh_points(in_mesh):
    """ Returns array mesh points from .vtk mesh

    Args:
        in_mesh (String): path to .vtk mesh containing points

    Returns:
        out_mesh_points (Numpy array): array of mesh points ordered in
            the same way as .vtk file points
    """
    import vtk
    import numpy as np
    # read .vtk file
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(in_mesh)
    reader.Update()
    polydata = reader.GetOutput()
    # populate numpy array
    point_number = polydata.GetNumberOfPoints()
    out_mesh_points = np.zeros((point_number, 3), np.float)
    for point_index in range(point_number):
        point_x, point_y, point_z = polydata.GetPoint(point_index)
        out_mesh_points[point_index, 0] = point_x
        out_mesh_points[point_index, 1] = point_y
        out_mesh_points[point_index, 2] = point_z

    return out_mesh_points



def generate_mesh_reoriented(
        in_mesh,
        in_mesh_reoriented_points,
        structure_name):
    """Reorient mesh with new coordinates

    Returns .vtk surface from 1) an original .vtk surface (to keep
    connectivity between vertices) and 2) an array of new point
    coordinates.

    Args:
        in_mesh (String): path to .vtk mesh containing the source
            surface. This is used to know the connectivity between
            vertices in the source mesh.
        in_mesh_reoriented_points (Numpy array): new point coordinates
            after rigid registration
        structure_name (String): name of the structure which is being
            saved. This is used to give a unique name to each structure
            mesh associated to the same subject.

    Returns:
        out_mesh_reoriented (String): path to .vtk mesh containing the
            rigidly registered surface
    """
    import os
    import vtk

    # read .vtk in_mesh
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(in_mesh)
    reader.Update()
    polydata = reader.GetOutput()
    # populate vtk points structure with the registered points
    vtk_reoriented_points = vtk.vtkPoints()
    point_number = polydata.GetNumberOfPoints()
    for point_index in range(point_number):
        reoriented_point = in_mesh_reoriented_points[point_index, :]
        vtk_reoriented_points.InsertNextPoint(
            reoriented_point[0],
            reoriented_point[1],
            reoriented_point[2])
    # add reoriented points to polydata
    polydata.SetPoints(vtk_reoriented_points)
    # write polydata to output file
    current_path = os.path.abspath('./')
    out_mesh_reoriented = '{0}/out_reoriented_{1}.vtk'.format(
        current_path,
        structure_name)
    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(polydata)
    writer.SetFileName(out_mesh_reoriented)
    writer.Update()

    return out_mesh_reoriented



class rigid_registration(object):
    """
    Entire class copied and pasted from Pycpd:
    'Pure Numpy Implementation of the Coherent Point Drift Algorithm'
    (https://github.com/siavashk/pycpd)
    """

    def __init__(self, X, Y, R=None, t=None, s=None, sigma2=None, maxIterations=100, tolerance=0.001, w=0):
        import numpy as np

        if X.shape[1] != Y.shape[1]:
          raise 'Both point clouds must have the same number of dimensions!'

        self.X             = X
        self.Y             = Y
        self.TY            = Y
        (self.N, self.D)   = self.X.shape
        (self.M, _)        = self.Y.shape
        self.R             = np.eye(self.D) if R is None else R
        self.t             = np.atleast_2d(np.zeros((1, self.D))) if t is None else t
        self.s             = 1 if s is None else s
        self.sigma2        = sigma2
        self.iteration     = 0
        self.maxIterations = maxIterations
        self.tolerance     = tolerance
        self.w             = w
        self.q             = 0
        self.err           = 0

    def register(self, callback):
        import numpy as np
        self.initialize()

        while self.iteration < self.maxIterations and self.err > self.tolerance:
            self.iterate()
            if callback:
                callback(iteration=self.iteration, error=self.err, X=self.X, Y=self.TY)

        return self.TY, self.R, np.atleast_2d(self.t), self.s

    def iterate(self):
        self.EStep()
        self.MStep()
        self.iteration += 1

    def MStep(self):
        self.updateTransform()
        self.transformPointCloud()
        self.updateVariance()

    def updateTransform(self):
        import numpy as np
        muX = np.divide(np.sum(np.dot(self.P, self.X), axis=0), self.Np)
        muY = np.divide(np.sum(np.dot(np.transpose(self.P), self.Y), axis=0), self.Np)

        self.XX = self.X - np.tile(muX, (self.N, 1))
        YY      = self.Y - np.tile(muY, (self.M, 1))

        self.A = np.dot(np.transpose(self.XX), np.transpose(self.P))
        self.A = np.dot(self.A, YY)

        U, _, V = np.linalg.svd(self.A, full_matrices=True)
        C = np.ones((self.D, ))
        C[self.D-1] = np.linalg.det(np.dot(U, V))

        self.R = np.dot(np.dot(U, np.diag(C)), V)

        self.YPY = np.dot(np.transpose(self.P1), np.sum(np.multiply(YY, YY), axis=1))

        self.s = np.trace(np.dot(np.transpose(self.A), self.R)) / self.YPY

        self.t = np.transpose(muX) - self.s * np.dot(self.R, np.transpose(muY))

    def transformPointCloud(self, Y=None):
        import numpy as np
        if Y is None:
            self.TY = self.s * np.dot(self.Y, np.transpose(self.R)) + np.tile(np.transpose(self.t), (self.M, 1))
            return
        else:
            return self.s * np.dot(Y, np.transpose(self.R)) + np.tile(np.transpose(self.t), (Y.shape[0], 1))

    def updateVariance(self):
        import numpy as np
        qprev = self.q

        trAR     = np.trace(np.dot(self.A, np.transpose(self.R)))
        xPx      = np.dot(np.transpose(self.Pt1), np.sum(np.multiply(self.XX, self.XX), axis =1))
        self.q   = (xPx - 2 * self.s * trAR + self.s * self.s * self.YPY) / (2 * self.sigma2) + self.D * self.Np/2 * np.log(self.sigma2)
        self.err = np.abs(self.q - qprev)

        self.sigma2 = (xPx - self.s * trAR) / (self.Np * self.D)

        if self.sigma2 <= 0:
            self.sigma2 = self.tolerance / 10

    def initialize(self):
        import numpy as np
        self.TY = self.s * np.dot(self.Y, np.transpose(self.R)) + np.repeat(self.t, self.M, axis=0)
        if not self.sigma2:
            XX = np.reshape(self.X, (1, self.N, self.D))
            YY = np.reshape(self.TY, (self.M, 1, self.D))
            XX = np.tile(XX, (self.M, 1, 1))
            YY = np.tile(YY, (1, self.N, 1))
            diff = XX - YY
            err  = np.multiply(diff, diff)
            self.sigma2 = np.sum(err) / (self.D * self.M * self.N)

        self.err  = self.tolerance + 1
        self.q    = -self.err - self.N * self.D/2 * np.log(self.sigma2)

    def EStep(self):
        import numpy as np
        P = np.zeros((self.M, self.N))

        for i in range(0, self.M):
            diff     = self.X - np.tile(self.TY[i, :], (self.N, 1))
            diff    = np.multiply(diff, diff)
            P[i, :] = P[i, :] + np.sum(diff, axis=1)

        c = (2 * np.pi * self.sigma2) ** (self.D / 2)
        c = c * self.w / (1 - self.w)
        c = c * self.M / self.N

        P = np.exp(-P / (2 * self.sigma2))
        den = np.sum(P, axis=0)
        den = np.tile(den, (self.M, 1))
        den[den==0] = np.finfo(float).eps

        self.P   = np.divide(P, den)
        self.Pt1 = np.sum(self.P, axis=0)
        self.P1  = np.sum(self.P, axis=1)
        self.Np  = np.sum(self.P1)



def rigid_register(source_points, target_points):
    """Rigid registration between a source and a target point set

    Args:
        source_points (Numpy array): [n1, 3] array
        target_points (Numpy array): [n2, 3] array

    Returns:
        registered_source_points (Numpy array): [n1, 3] transformation
            matrix
    """
    from functools import partial
    # check input is correct
    if source_points.shape[1] != 3:
        raise RuntimeError("Wrong dimensions for source_points")
    if target_points.shape[1] != 3:
        raise RuntimeError("Wrong dimensions for target_points")

    # define empty visualise function
    def visualise(iteration, error, X, Y, ax):
        """ Function used to create an 'empty' visualiser that will
        be passed to the Pycpd rigid_registration object
        """
        return

    # create rigid registration object
    max_iterations = 100
    error_tolerance = 0.001
    reg = rigid_registration(
        target_points,
        source_points,
        maxIterations=max_iterations,
        tolerance=error_tolerance)
    call_back = partial(visualise, ax=None)
    reg.register(call_back)
    registered_source_points = reg.TY

    return registered_source_points



def rigid_register_colin27(
        in_persubject_mesh_list,
        in_colin27_mesh,
        in_structure_list):
    """Register structures to corresponding Colin27 structures

    For a given subject, register a set of structures (combined as a
    group) to the corresponding set of Colin27 structures. Apply the
    rigid registration to each structures in the set and output the
    registered structures.
    This re-uses the code from Pycpd: 'Pure Numpy Implementation of the
    Coherent Point Drift Algorithm' (https://github.com/siavashk/pycpd)

    Args:
        in_persubject_mesh_list (list of String): list of per-subject
            .vtk mesh paths (as many surfaces per subjects as there are
            structures)
        in_colin27_mesh (String): path to colin 27 .vtk corresponding
            structure
        in_structure_list (list of String): names of the structures that
            are being processed

    Returns:
        out_persubject_regmesh_list (list of String): list of
            per-subject .vtk rigidly registered meshs
    """
    import numpy as np
    import deformetrica_preprocessing_utils as utils

    out_persubject_regmesh_list = []
    subject_overall_points = np.zeros((0, 3), np.float)
    structure_boundaries_list = []
    colin27_overall_points = np.zeros((0, 3), np.float)
    structure_number = len(in_persubject_mesh_list)
    # assemble subject (source) and Colin 27 (target) point sets
    for structure_index in range(structure_number):
        # read meshes
        # subject structure
        subject_structure_mesh = in_persubject_mesh_list[structure_index]
        subject_structure_points = utils.get_mesh_points(subject_structure_mesh)
        # Colin 27 structure
        colin27_structure_mesh = in_colin27_mesh[structure_index]
        colin27_structure_points = utils.get_mesh_points(colin27_structure_mesh)
        # store boundaries of the current structure block within the overall set
        block_start = subject_overall_points.shape[0]
        block_length = subject_structure_points.shape[0]
        block_end = block_start + block_length
        structure_boundaries_list.append([block_start, block_end])
        # merge meshes
        # subject
        subject_overall_points = np.concatenate(
            (subject_overall_points, subject_structure_points), axis=0)
        # Colin27
        colin27_overall_points = np.concatenate(
            (colin27_overall_points, colin27_structure_points), axis=0)
    # rigidly register subject to Colin27 point sets
    rigid_subject_overall_points = utils.rigid_register(
        subject_overall_points,
        colin27_overall_points)
    # generate reoriented surfaces for each of the structures associated to the
    # subject
    for structure_index in range(structure_number):
        # read subject current structure mesh
        subject_structure_mesh = in_persubject_mesh_list[structure_index]
        # extract rigidly registered points for current structure
        structure_boundaries = structure_boundaries_list[structure_index]
        block_start = structure_boundaries[0]
        block_end = structure_boundaries[1]
        rigid_subject_structure_points = rigid_subject_overall_points[
            block_start:block_end]
        # generate mesh for current structure
        structure_name = in_structure_list[structure_index]
        out_rigid_subject_structure_mesh = utils.generate_mesh_reoriented(
            subject_structure_mesh,
            rigid_subject_structure_points,
            structure_name)
        out_persubject_regmesh_list.append(out_rigid_subject_structure_mesh)

    return out_persubject_regmesh_list



def link_basedirs_structures(in_basedir_list, in_structure_list):
    """Links list of datasink base directories to list of structures IDs

    Args:
        in_basedir_list (list of String): list of datasink base
            directory paths
        in_structure_list (list of String): list of structures

    Returns:
        out_basedir_list (list of String): reordered list of datasink
            base directory paths
    """
    import deformetrica_preprocessing_utils as utils
    [out_basedir_list, dummy, dummy] = utils.link_objects1_objects2(
        in_basedir_list,
        in_structure_list)

    return out_basedir_list



def initialise_mapnode_datasink_regexp(datasink_name):
    """Create regexp to get rid of unwanted mapnode datasink folders

    All mapnode datasinks seem to automatically create the same type of
    (undesirable) folders. This function will produce the regular
    expressions to eliminate those.
    Eliminate 'trait_added' and [datasink_name] folders.

    Args:
        datasink_name (String): name of the datasink

    Returns:
        datasink_regexp_substitution_list (list of tuples): list of
            regular expressions to remove folders that have been
            created by a mapnode datasink
    """
    datasink_regexp_substitution_list = [
        (r'trait_added', r''),
        (r'{0}.*'.format(datasink_name), r'')
        ]

    return datasink_regexp_substitution_list
