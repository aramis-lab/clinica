"""

"""

def get_class_images(class_images, index_list):
    """
    utility method to extract class images
    from a multi session class_images set:
    class nb : tissue type
    1 : Grey Matter
    2 : White Matter
    3 : CerebroSpinal Fluid
    4 : Skull
    5 : Out-of-brain soft tissue
    6 : Head surrounding

    :param class_images: image set from which to extract images.
    :param index_list: index list of the classes to extract.
    :return: extracted images in a list of lists (without empty lists).

    Example
    -------
    >>> class_n_images = get_class_images(class_images,[1,2])
    """

    # Declare class images list
    class_n_images = {}
    for idx in index_list:
        class_n_images[idx] = []

    for session in class_images:
        for idx in index_list:
            class_n_images[idx].extend(session[idx-1])

    result = []
    for idx in index_list:
        if class_n_images[idx]:
            result.append(class_n_images[idx])

    return result


def step2(in_hello_word):
    """Example function for Step 2.
    """

    print(in_hello_word + " from the step 2 of the T1 SPM Dartel Clinica pipeline.")