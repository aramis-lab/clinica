from typing import Callable, Dict


class Query:
    default_queries = {}
    repr_indent = 4

    def __init__(self, query: Dict):
        self.query = self.format_query(query)

    def format_query(self, query: Dict):
        formatted_query = {}
        for k, q in query.items():
            if k in self.default_queries:
                formatted_query[k] = self.combine_queries(self.default_queries[k], q)
        return formatted_query

    @staticmethod
    def combine_queries(default_query: Dict, user_query: Dict):
        pass

    def __len__(self):
        return len(self.query)

    def __str__(self):
        import json

        return json.dumps(self.query, indent=self.repr_indent)


class BIDSQuery(Query):
    default_queries = {
        "T1w": {"datatype": "anat", "suffix": "T1w", "extension": [".nii.gz"]}
    }

    @staticmethod
    def combine_queries(default_query: Dict, user_query: Dict) -> Dict:
        return {**default_query, **user_query}


class CAPSQuery(Query):
    @staticmethod
    def combine_queries(query_maker: Callable, user_query: Dict) -> Dict:
        return query_maker(**user_query)


class CAPSFileQuery(CAPSQuery):
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
    from clinica.utils.input_files import t1_volume_final_group_template

    default_queries = {
        "dartel_template": t1_volume_final_group_template,
    }


def query_factory(query: dict, query_type: str) -> Query:
    from clinica.pydra.query import BIDSQuery, CAPSFileQuery, CAPSGroupQuery

    if query_type == "bids":
        return BIDSQuery(query)
    elif query_type == "caps_file":
        return CAPSFileQuery(query)
    elif query_type == "caps_group":
        return CAPSGroupQuery(query)
    else:
        raise ValueError(".")
