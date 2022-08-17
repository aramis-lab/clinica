def test_ApplySegmentationDeformation():
    """Uses SPM to apply a deformation field obtained from Segmentation routine to a given file"""

    import clinica.pydra.t1_volume.tasks as t1vol_tasks

    inv = t1vol_tasks.ApplySegmentationDeformation()
    assert inv._jobname == "defs"
    assert inv._jobtype == "util"
    assert inv.input_spec().get_traitsfree() == {
        "mask": 0,
        "mfile": True,
        "use_v8struct": True,
    }
    return
