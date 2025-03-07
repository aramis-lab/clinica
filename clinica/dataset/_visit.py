from dataclasses import dataclass

__all__ = ["Visit"]


@dataclass(frozen=True)
class Visit:
    subject: str
    session: str

    def __lt__(self, obj):
        return (self.subject < obj.subject) or (
            self.subject == obj.subject and self.session < obj.session
        )

    def __gt__(self, obj):
        return (self.subject > obj.subject) or (
            self.subject == obj.subject and self.session > obj.session
        )

    def __str__(self) -> str:
        return f"{self.subject} {self.session}"
