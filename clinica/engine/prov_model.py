from attr import define, field
import attr
import typing
from typing import Union, List
from abc import ABC, abstractmethod


#  Define PROV abstract concepts


@define
class Identifier:
    id: int


class ProvElement(ABC):
    @property
    @classmethod
    @abstractmethod
    def id(cls):
        """id is required for ProvElements"""
        return NotImplementedError

    @property
    def attributes(cls):
        """attributes are optional"""
        return NotImplementedError


class ProvRelation(ABC):

    id: Identifier
    src: ProvElement
    dest: ProvElement


# Define PROV Types


@define
class ProvEntity(ProvElement):
    """Provenance Entity element"""

    id: Identifier = field(validator=[attr.validators.instance_of(Identifier)])
    attributes: dict


@define
class ProvActivity(ProvElement):
    """Provenance Activity element"""

    id: Identifier = field(validator=[attr.validators.instance_of(Identifier)])
    attributes: dict


@define
class ProvAgent(ProvElement):
    """Provenance Agent element"""

    id: Identifier = field(validator=[attr.validators.instance_of(Identifier)])
    attributes: dict


# Define PROV Relations


@define
class ProvGeneration(ProvRelation):
    id: Identifier = field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of(Identifier)),
    )

    src: ProvActivity = field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of(ProvActivity)),
    )
    dest: ProvEntity = field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of(ProvEntity)),
    )

    # entity: an identifier (e) for a created entity;
    # activity: an OPTIONAL identifier (a) for the activity that creates the entity;
    # time: an OPTIONAL "generation time" (t), the time at which the entity was completely created;
    # attributes: an OPTIONALa


@define
class ProvUsage(ProvRelation):
    pass


@define
class ProvAssociation(ProvRelation):
    pass


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

    context: dict
    entries: list[ProvEntry]
