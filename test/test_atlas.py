from nose.tools import assert_equals
from nose.tools import assert_true

from clinica.utils.atlas import *

def test_AtlasName():
    assert_equals(Atlas_JHUDTI81_1mm().get_name_atlas(), "JHUDTI81")
    assert_equals(Atlas_JHUDTI81_2mm().get_name_atlas(), "JHUDTI81")
    assert_equals(Atlas_JHUTracts0_1mm().get_name_atlas(), "JHUTracts0")
    assert_equals(Atlas_JHUTracts0_2mm().get_name_atlas(), "JHUTracts0")
    assert_equals(Atlas_JHUTracts25_1mm().get_name_atlas(), "JHUTracts25")
    assert_equals(Atlas_JHUTracts25_2mm().get_name_atlas(), "JHUTracts25")
    assert_equals(Atlas_JHUTracts50_1mm().get_name_atlas(), "JHUTracts50")
    assert_equals(Atlas_JHUTracts50_2mm().get_name_atlas(), "JHUTracts50")

def test_AtlasResolution():
    assert_equals(Atlas_JHUDTI81_1mm().get_spatial_resolution(), "1x1x1")
    assert_equals(Atlas_JHUDTI81_2mm().get_spatial_resolution(), "2x2x2")
    assert_equals(Atlas_JHUTracts0_1mm().get_spatial_resolution(), "1x1x1")
    assert_equals(Atlas_JHUTracts0_2mm().get_spatial_resolution(), "2x2x2")
    assert_equals(Atlas_JHUTracts25_1mm().get_spatial_resolution(), "1x1x1")
    assert_equals(Atlas_JHUTracts25_2mm().get_spatial_resolution(), "2x2x2")
    assert_equals(Atlas_JHUTracts50_1mm().get_spatial_resolution(), "1x1x1")
    assert_equals(Atlas_JHUTracts50_2mm().get_spatial_resolution(), "2x2x2")

def test_AtlasLabels():
    import os
    assert_true(os.path.isfile(Atlas_JHUDTI81_1mm().get_atlas_labels()))
    assert_true(os.path.isfile(Atlas_JHUDTI81_2mm().get_atlas_labels()))
    assert_true(os.path.isfile(Atlas_JHUTracts0_1mm().get_atlas_labels()))
    assert_true(os.path.isfile(Atlas_JHUTracts0_2mm().get_atlas_labels()))
    assert_true(os.path.isfile(Atlas_JHUTracts25_1mm().get_atlas_labels()))
    assert_true(os.path.isfile(Atlas_JHUTracts25_2mm().get_atlas_labels()))
    assert_true(os.path.isfile(Atlas_JHUTracts50_1mm().get_atlas_labels()))
    assert_true(os.path.isfile(Atlas_JHUTracts50_2mm().get_atlas_labels()))

def test_AtlasMap():
    import os
    assert_true(os.path.isfile(Atlas_JHUDTI81_1mm().get_atlas_map()))
    assert_true(os.path.isfile(Atlas_JHUDTI81_2mm().get_atlas_map()))
    assert_true(os.path.isfile(Atlas_JHUTracts0_1mm().get_atlas_map()))
    assert_true(os.path.isfile(Atlas_JHUTracts0_2mm().get_atlas_map()))
    assert_true(os.path.isfile(Atlas_JHUTracts25_1mm().get_atlas_map()))
    assert_true(os.path.isfile(Atlas_JHUTracts25_2mm().get_atlas_map()))
    assert_true(os.path.isfile(Atlas_JHUTracts50_1mm().get_atlas_map()))
    assert_true(os.path.isfile(Atlas_JHUTracts50_2mm().get_atlas_map()))

