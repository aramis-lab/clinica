import abc
from typing import Callable, Dict, Optional


class Query:
    """Base Query class.

    Attributes
    ----------
    query : Dict, optional
        Raw input query defined as a dictionary.
        Default=None.
    """

    default_queries = {}
    repr_indent = 4

    def __init__(self, query: Optional[Dict] = None):
        self.query = self.format_query(query)

    def format_query(self, query: Optional[Dict] = None) -> Dict:
        """Format the input query by combining and filtering with
        the class default queries available.

        Parameters
        ----------
        query : Dict
            Raw input query to format.

        Returns
        -------
        formatted_query : Dict
            Resulting formatted query.
        """
        formatted_query = {}
        if not query:
            return formatted_query
        for k, q in query.items():
            if k in self.default_queries:
                formatted_query[k] = self.combine_queries(self.default_queries[k], q)
        return formatted_query

    @staticmethod
    @abc.abstractmethod
    def combine_queries(default_query: Dict, user_query: Dict) -> Dict:
        pass

    def __len__(self):
        return len(self.query)

    def __str__(self):
        import json

        return json.dumps(self.query, indent=self.repr_indent)


class BIDSQuery(Query):
    """BIDSQuery class."""

    default_queries = {
        "T1w": {"datatype": "anat", "suffix": "T1w", "extension": [".nii.gz"]}
    }

    @staticmethod
    def combine_queries(default_query: Dict, user_query: Dict) -> Dict:
        """Method combining a default query and a user-provided-query.

        Parameters
        -----------
        default_query : Dict
            The default query to use for combining.

        user_query : Dict
            The user query to use for combining.

        Returns
        -------
        Dict :
            The combined query.
        """
        return {**default_query, **user_query}


class CAPSQuery(Query):
    """CAPSQuery class."""

    @staticmethod
    def combine_queries(query_maker: Callable, user_query: Dict) -> Dict:
        """Method building a proper CAPS query from a query maker
        function and a user-provided-query.

        Parameters
        -----------
        query_maker : Callable
            The function to use for building the CAPS query.
            These functions are already implemented in Clinica.
            They just need to be imported and linked to the proper key.

        user_query : Dict
            The user query to use. This is basically a kwargs dictionary
            to pass to the `query_maker` in order to build the query with
            the correct inputs.

        Returns
        -------
        Dict :
            The combined query.
        """
        return query_maker(**user_query)


class CAPSFileQuery(CAPSQuery):
    """CAPSFileQuery class."""

    from clinica.utils.input_files import (
        t1_volume_deformation_to_template,
        t1_volume_native_tpm,
        t1_volume_native_tpm_in_mni,
    )

    default_queries = {
        "mask_tissues": t1_volume_native_tpm_in_mni,
        "flow_fields": t1_volume_deformation_to_template,
        "pvc_mask_tissues": t1_volume_native_tpm,
    }


class CAPSGroupQuery(CAPSQuery):
    """CAPSGroupQuery class."""

    from clinica.utils.input_files import t1_volume_final_group_template

    default_queries = {
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
