import abc
from typing import Callable, Dict, Optional


class Query:
    """Base Query class.

    Examples
    --------
    >>> from clinica.pydra.query import Query
    >>> q = Query()
    >>> len(q)
    0
    """

    _default_queries = {}
    repr_indent = 4

    def __init__(self, query: Optional[Dict] = None):
        self._query = self.format_query(query)

    @property
    def query(self):
        return self._query

    def format_query(self, input_query: Optional[Dict] = None) -> Dict:
        """Format the input query in a form suitable for DataGrabbers.

        If the query is None, return an empty dictionary.

        Otherwise, `format_query` does two things:

            - parse the input query
            - potentially combine the parsed input query with default values

        The input query is a dictionary but has a heterogeneous structure
        depending on whether we wish to query BIDS or CAPS folders.

        More precisely, it can have the following form :

            {"label" : {query dict, possibly incomplete}}

            example : {"T1w": {"suffix": "T1w"}}

        which is the usual format for BIDS queries. In this case, we simply need
        to filter out queries for labels not supported for BIDS (for example we'd
        keep "T1w" but drop "dartel" or "mask_tissues" if these were present since
        these are labels for CAPS queries).

        The input query can also have the following form:

            {"label" : {kwargs for a query maker}}

            example: {"mask_tissues": {'tissue_number': (1, 2), 'modulation': False}}

        which is the usual format for CAPS queries. In this case, we need to
        filter out queries for labels not supported (for example we'd keep "mask_tissues"
        but drop "T1w" or "pet" if these were present since these are labels for BIDS queries).
        And we need to call the appropriate query maker to obtain a proper query dict.
        In this example, the query maker corresponding to the label "mask_tissues" is the
        function `clinica.utils.input_files.t1_volume_native_tpm_in_mni`. We call this
        function with arguments `tissue_number=(1, 2), modulation=False`.

        After parsing, the query structure is the same for both BIDS and CAPS:

            {"label" : {query dict, possibly incomplete} or list of query dicts}

        The query is then completed with potential default values:

            {"label" : {query dict} or list of query dicts}

            example: {"T1w": {"suffix": "T1w", "extension": [".nii.gz"]}}

        Parameters
        ----------
        input_query : Dict, optional
            Raw input query to format.

        Returns
        -------
        formatted_query : Dict
            Resulting formatted query.
        """
        formatted_query = {}
        if not input_query:
            return formatted_query
        for k, q in self.parse_query(input_query).items():
            if isinstance(q, dict):
                formatted_query[k] = {**self.default_query(k), **q}
            elif isinstance(q, list):
                formatted_query[k] = [{**self.default_query(k), **qq} for qq in q]
            else:
                raise TypeError(
                    f"Unexpected type {type(q)} for query {q}."
                    "Only dictionaries and lists are supported."
                )
        return formatted_query

    @abc.abstractmethod
    def parse_query(self, query: Dict) -> Dict:
        pass

    def default_query(self, label: str) -> Dict:
        """Returns the default query corresponding to the given label.
        If there is no such query, returns an empty dict.

        Parameters
        ----------
        label : str
            Label of the default query to look for.

        Returns
        -------
        Dict :
            The default query corresponding to the provided label.
        """
        return self._default_queries.get(label, {})

    def __getitem__(self, name):
        return self._query[name]

    def __iter__(self):
        return iter(self._query)

    def keys(self):
        return self._query.keys()

    def items(self):
        return self._query.items()

    def values(self):
        return self._query.values()

    def __len__(self):
        return len(self._query)

    def __str__(self):
        import json

        return json.dumps(self._query, indent=self.repr_indent)


class BIDSQuery(Query):
    """BIDSQuery class.

    Examples
    --------
    >>> from clinica.pydra.query import BIDSQuery
    >>> q = BIDSQuery({"T1w": {"suffix": "T1w"}})
    >>> len(q)
    1
    >>> q.query
    {'T1w': {'suffix': 'T1w', 'extension': ['.nii.gz']}}
    """

    _default_queries = {
        "T1w": {"suffix": "T1w", "extension": [".nii.gz"]},
        "pet": {"suffix": "pet", "extension": [".nii.gz"]},
    }

    def parse_query(self, query: Dict) -> Dict:
        """Parse the provided query.

        BIDS query are already specified in a proper dictionary format.
        Parsing is therefore equivalent to filtering out queries corresponding
        to invalid labels (for example labels for CAPS).
        """
        return {k: v for k, v in query.items() if k in self._default_queries}


class CAPSQuery(Query):
    """CAPSQuery class."""

    # Query makers are specified in subclasses
    _query_makers = {}

    def parse_query(self, query: Dict) -> Dict:
        """Parse the provided query.

        For CAPS, the user is providing a dictionary of kwargs for some function
        responsible for building the query dictionary (for example
        `clinica.utils.input_files.t1_volume_deformation_to_template`).
        These functions are stored in a mapping and accessed through the `query_maker`
        method which uses the labels to look for these functions.
        Parsing the query is equivalent to calling these functions for each label with
        the associated user-provided kwargs.

        Parameters
        ----------
        query : Dict
            The user provided query to be parsed.

        Returns
        -------
        parsed_query : Dict
            The parsed query.
        """
        parsed_query = {}
        for label, params in query.items():
            query_maker = self._query_maker(label)
            formatted_query = query_maker(**params)
            if len(formatted_query) > 0:
                parsed_query[label] = formatted_query
        return parsed_query

    def _query_maker(self, label: str) -> Callable:
        """Retrieve the query-making-function associated with the
        provided label.

        Parameters
        ----------
        label : str
            The label to use to get the function.

        Returns
        -------
        Callable :
            The function to use to generate the query dictionary for this label.
            If the label does not match any entry, a default maker which return
             an empty dict for any passed parameters is returned.
        """
        return self._query_makers.get(label, lambda **kwargs: {})


class CAPSFileQuery(CAPSQuery):
    """CAPSFileQuery class.

    This class only holds the mapping between labels and query makers
    to be used with the `CAPSFileDataGrabber`.

    Examples
    --------
    >>> from clinica.pydra.query import CAPSFileQuery
    >>> q = CAPSFileQuery({'mask_tissues': {'tissue_number': (1, 2), 'modulation': False}})
    >>> len(q)
    1
    >>> import json
    >>> print(json.dumps(q.query, indent=2))
    {
        "mask_tissues": [
            {
                "pattern": "t1/spm/segmentation/normalized_space/*_*_T1w_segm-graymatter_space-Ixi549Space_modulated-off_probability.nii*",
                "description": "Tissue probability map graymatter based on native MRI in MNI space (Ixi549) without modulation.",
                "needed_pipeline": "t1-volume-tissue-segmentation"
            },
            {
                "pattern": "t1/spm/segmentation/normalized_space/*_*_T1w_segm-whitematter_space-Ixi549Space_modulated-off_probability.nii*",
                "description": "Tissue probability map whitematter based on native MRI in MNI space (Ixi549) without modulation.",
                "needed_pipeline": "t1-volume-tissue-segmentation"
            }
        ]
    }
    """

    from clinica.utils.input_files import (
        t1_volume_dartel_input_tissue,
        t1_volume_deformation_to_template,
        t1_volume_native_tpm,
        t1_volume_native_tpm_in_mni,
    )

    _query_makers = {
        "mask_tissues": t1_volume_native_tpm_in_mni,
        "flow_fields": t1_volume_deformation_to_template,
        "pvc_mask_tissues": t1_volume_native_tpm,
        "dartel_input_tissue": t1_volume_dartel_input_tissue,
    }


class CAPSGroupQuery(CAPSQuery):
    """CAPSGroupQuery class.

    This class only holds the mapping between the labels and query makers
    to be used with the `CAPSGroupDataGrabber`.

    Examples
    --------
    >>> from clinica.pydra.query import CAPSGroupQuery
    >>> import json
    >>> q = CAPSGroupQuery({'dartel_template': {'group_label': "UnitTest"}})
    >>> len(q)
    1
    >>> print(json.dumps(q.query, indent=2))
    {
        "dartel_template": {
        "pattern": "group-UnitTest/t1/group-UnitTest_template.nii*",
        "description": "T1w template file of group UnitTest",
        "needed_pipeline": "t1-volume or t1-volume-create-dartel"
        }
    }
    """

    from clinica.utils.input_files import t1_volume_final_group_template

    _query_makers = {
        "dartel_template": t1_volume_final_group_template,
    }


def query_factory(query: Dict, query_type: str) -> Query:
    """Return the Query instance initialized from passed `query`,
    based on provided `query_type`.

    Parameters
    ----------
    query : dict
        The query in dictionary format.

    query_type : {"bids", "caps_file", "caps_group"}
        The type of query to build.

    Returns
    -------
    Query :
        The resulting Query instance.

    Raises
    ------
    ValueError
        If `query_type` is not a valid value.
    """
    if query_type == "bids":
        return BIDSQuery(query)
    elif query_type == "caps_file":
        return CAPSFileQuery(query)
    elif query_type == "caps_group":
        return CAPSGroupQuery(query)
    else:
        raise ValueError(
            f"query_factory() received an invalid argument for 'query_type' : {query_type}"
            f"Supported values are 'bids', 'caps_file', and 'caps_group'."
        )
