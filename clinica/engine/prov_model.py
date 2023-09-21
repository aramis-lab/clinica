from abc import ABC, abstractmethod
from typing import List

import attr
from attr import define, field

#  Define PROV abstract concepts


@define
class ProvContext:
    _namespaces: list


@define
class Namespace:
    id: str
    uri: str


@define
class Identifier:
    label: str = field(
        validator=attr.validators.optional(attr.validators.instance_of(str)),
    )

    def __repr__(self):
        return "%s" % self.label


class ProvElement(ABC):
    @property
    @classmethod
    @abstractmethod
    def uid(cls):
        """id is required for ProvElements"""
        return NotImplementedError

    @property
    def attributes(cls):
        """attributes are optional"""
        return NotImplementedError

    @classmethod
    def get_type(cls):
        return type(cls).__name__


class ProvRelation(ABC):

    id: Identifier
    src: ProvElement
    dest: ProvElement


# Define PROV Types


@define
class ProvEntity(ProvElement):
    """Provenance Entity element"""

    uid: Identifier = field(validator=[attr.validators.instance_of(Identifier)])
    attributes: dict = field(default=attr.Factory(dict))

    def unstruct(self):
        return {"id": str(self.uid), **self.attributes}


@define
class ProvActivity(ProvElement):
    """Provenance Activity element"""

    uid: Identifier = field(validator=[attr.validators.instance_of(Identifier)])
    attributes: dict = field(default=attr.Factory(dict))

    def unstruct(self):
        return {"id": str(self.uid), **self.attributes}


@define
class ProvAgent(ProvElement):
    """Provenance Agent element"""

    uid: Identifier = field(validator=[attr.validators.instance_of(Identifier)])
    attributes: dict = field(default=attr.Factory(dict))

    def unstruct(self):
        return {"id": str(self.uid), **self.attributes}


@define
class ProvEntry:
    """
    A prov entry in triple form
    """

    subject: ProvElement
    predicate: ProvRelation
    object: ProvElement


@define
class ProvRecord:
    """
    A provenance document containting a PROV context and a list of entries
    """

    context: ProvContext = field()
    elements: List[ProvElement] = field(default=[])

    def __getitem__(self, idx):
        for element in self.elements:
            if element.uid == idx:
                return element

    def json(self):

        json_dict = {}

        context_keys = [x.id for x in self.context._namespaces]
        context_vals = [y.uri for y in self.context._namespaces]

        json_dict["@context"] = dict(zip(context_keys, context_vals))

        json_dict["prov:Agent"] = [
            x.unstruct() for x in self.elements if isinstance(x, ProvAgent)
        ]
        json_dict["prov:Activity"] = [
            x.unstruct() for x in self.elements if isinstance(x, ProvActivity)
        ]
        json_dict["prov:Entity"] = [
            x.unstruct() for x in self.elements if isinstance(x, ProvEntity)
        ]
        return json_dict
