from nose.tools import assert_equals
from nose.tools import assert_true

from clinica.utils.atlas import *

def test_AtlasName():
    assert_equals(JHUDTI81_1mm().get_name_atlas(), "JHUDTI81")
    assert_equals(JHUDTI81_2mm().get_name_atlas(), "JHUDTI81")
    assert_equals(JHUTracts0_1mm().get_name_atlas(), "JHUTracts0")
    assert_equals(JHUTracts0_2mm().get_name_atlas(), "JHUTracts0")
    assert_equals(JHUTracts25_1mm().get_name_atlas(), "JHUTracts25")
    assert_equals(JHUTracts25_2mm().get_name_atlas(), "JHUTracts25")
    assert_equals(JHUTracts50_1mm().get_name_atlas(), "JHUTracts50")
    assert_equals(JHUTracts50_2mm().get_name_atlas(), "JHUTracts50")
    assert_equals(AAL2().get_name_atlas(), "AAL2")
    assert_equals(Hammers().get_name_atlas(), "Hammers")
    assert_equals(Neuromorphometrics().get_name_atlas(), "Neuromorphometrics")
    assert_equals(AICHA().get_name_atlas(), "AICHA")
    assert_equals(LPBA40().get_name_atlas(), "LPBA40")


def test_AtlasResolution():
    assert_equals(JHUDTI81_1mm().get_spatial_resolution(), "1x1x1")
    assert_equals(JHUDTI81_2mm().get_spatial_resolution(), "2x2x2")
    assert_equals(JHUTracts0_1mm().get_spatial_resolution(), "1x1x1")
    assert_equals(JHUTracts0_2mm().get_spatial_resolution(), "2x2x2")
    assert_equals(JHUTracts25_1mm().get_spatial_resolution(), "1x1x1")
    assert_equals(JHUTracts25_2mm().get_spatial_resolution(), "2x2x2")
    assert_equals(JHUTracts50_1mm().get_spatial_resolution(), "1x1x1")
    assert_equals(JHUTracts50_2mm().get_spatial_resolution(), "2x2x2")
    assert_equals(AAL2().get_spatial_resolution(), "1.5x1.5x1.5")
    assert_equals(Hammers().get_spatial_resolution(), "1.5x1.5x1.5")
    assert_equals(Neuromorphometrics().get_spatial_resolution(), "1.5x1.5x1.5")
    assert_equals(LPBA40().get_spatial_resolution(), "1.5x1.5x1.5")
    assert_equals(AICHA().get_spatial_resolution(), "1.5x1.5x1.5")




def test_AtlasLabels():
    import os
    assert_true(os.path.isfile(JHUDTI81_1mm().get_atlas_labels()))
    assert_true(os.path.isfile(JHUDTI81_2mm().get_atlas_labels()))
    assert_true(os.path.isfile(JHUTracts0_1mm().get_atlas_labels()))
    assert_true(os.path.isfile(JHUTracts0_2mm().get_atlas_labels()))
    assert_true(os.path.isfile(JHUTracts25_1mm().get_atlas_labels()))
    assert_true(os.path.isfile(JHUTracts25_2mm().get_atlas_labels()))
    assert_true(os.path.isfile(JHUTracts50_1mm().get_atlas_labels()))
    assert_true(os.path.isfile(JHUTracts50_2mm().get_atlas_labels()))
    assert_true(os.path.isfile(AAL2().get_atlas_labels()))
    assert_true(os.path.isfile(Hammers().get_atlas_labels()))
    assert_true(os.path.isfile(Neuromorphometrics().get_atlas_labels()))
    assert_true(os.path.isfile(AICHA().get_atlas_labels()))
    assert_true(os.path.isfile(LPBA40().get_atlas_labels()))




def test_AtlasMap():
    import os
    assert_true(os.path.isfile(JHUDTI81_1mm().get_atlas_map()))
    assert_true(os.path.isfile(JHUDTI81_2mm().get_atlas_map()))
    assert_true(os.path.isfile(JHUTracts0_1mm().get_atlas_map()))
    assert_true(os.path.isfile(JHUTracts0_2mm().get_atlas_map()))
    assert_true(os.path.isfile(JHUTracts25_1mm().get_atlas_map()))
    assert_true(os.path.isfile(JHUTracts25_2mm().get_atlas_map()))
    assert_true(os.path.isfile(JHUTracts50_1mm().get_atlas_map()))
    assert_true(os.path.isfile(JHUTracts50_2mm().get_atlas_map()))
    assert_true(os.path.isfile(AAL2().get_atlas_map()))
    assert_true(os.path.isfile(Neuromorphometrics().get_atlas_map()))
    assert_true(os.path.isfile(Hammers().get_atlas_map()))
    assert_true(os.path.isfile(AICHA().get_atlas_map()))
    assert_true(os.path.isfile(LPBA40().get_atlas_map()))


def test_Index():
    assert_equals(str(len(AAL2().get_index())), '121')
    assert_equals(str(len(Neuromorphometrics().get_index())), '141')
    assert_equals(str(len(Hammers().get_index())), '69')
    assert_equals(str(len(LPBA40().get_index())), '57')
    assert_equals(str(len(AICHA().get_index())), '385')


def test_AtlasTSVROI():
    import os
    assert_true(os.path.isfile(JHUDTI81_1mm().get_tsv_roi()))
    assert_true(os.path.isfile(JHUDTI81_2mm().get_tsv_roi()))
    assert_true(os.path.isfile(JHUTracts0_1mm().get_tsv_roi()))
    assert_true(os.path.isfile(JHUTracts0_2mm().get_tsv_roi()))
    assert_true(os.path.isfile(JHUTracts25_1mm().get_tsv_roi()))
    assert_true(os.path.isfile(JHUTracts25_2mm().get_tsv_roi()))
    assert_true(os.path.isfile(JHUTracts50_1mm().get_tsv_roi()))
    assert_true(os.path.isfile(JHUTracts50_2mm().get_tsv_roi()))
    assert_true(os.path.isfile(AAL2().get_tsv_roi()))
    assert_true(os.path.isfile(AICHA().get_tsv_roi()))
    assert_true(os.path.isfile(Hammers().get_tsv_roi()))
    assert_true(os.path.isfile(LPBA40().get_tsv_roi()))
    assert_true(os.path.isfile(Neuromorphometrics().get_tsv_roi()))

