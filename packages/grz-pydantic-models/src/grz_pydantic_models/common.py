from pydantic import BaseModel, ConfigDict
from pydantic.alias_generators import to_camel


class StrictBaseModel(BaseModel):
    model_config = ConfigDict(
        extra="forbid",
        validate_assignment=True,
        use_enum_values=True,
        alias_generator=to_camel,
    )
