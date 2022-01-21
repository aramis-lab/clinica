from attr import define, field
import attr
import typing
from typing import Union, List
from abc import ABC, abstractmethod


#  Define PROV abstract concepts


@define
class ProvContext:
    label: str
    link: str


@define
class Identifier:
    seed: int = field()
    label: str = field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of(str)),
    )


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

    id: Identifier = field(
        init=False, validator=[attr.validators.instance_of(Identifier)]
    )
    attributes: dict = field(default={})

    def __attrs_post_init__(self):
        self.id = Identifier(seed=0)


@define
class ProvActivity(ProvElement):
    """Provenance Activity element"""

    id: Identifier = field(
        init=False, validator=[attr.validators.instance_of(Identifier)]
    )
    attributes: dict = field(default={})

    def __attrs_post_init__(self):
        self.id = Identifier(seed=0)


@define
class ProvAgent(ProvElement):
    """Provenance Agent element"""

    id: Identifier = field(
        init=False, validator=[attr.validators.instance_of(Identifier)]
    )
    attributes: dict = field(
        default={},
        validator=attr.validators.optional(attr.validators.instance_of(dict)),
    )

    def __attrs_post_init__(self):
        self.id = Identifier(seed=0)


# Define PROV Relations


@define
class ProvGeneration(ProvRelation):
    id: Identifier = field(
        init=False,
        validator=attr.validators.optional(attr.validators.instance_of(Identifier)),
    )

    src: ProvActivity = field(
        init=False,
        validator=attr.validators.optional(attr.validators.instance_of(ProvActivity)),
    )
    dest: ProvEntity = field(
        init=False,
        validator=attr.validators.optional(attr.validators.instance_of(ProvEntity)),
    )

    def __attrs_post_init__(self):
        self.id = Identifier(seed=0)
        self.src = ProvActivity()
        self.dest = ProvEntity()

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

    context: dict = field(default={})
    entries: List[ProvEntry] = field(default=[])
