# Developer guide to migrations

## Create a new migration

From `packages/grz-db`:

```
uv run alembic revision -m "description of migration"
```

A new migration script is placed in `versions/`.
All you have to do is implement `upgrade()` to change the database as needed for the migration, including adding necessary `import`s at the top of the script as needed.
The migration metadata fields, such as `down_revision`, are populated by Alembic and shouldn't be changed.
The available operations can be browsed in the [Alembic documentation](https://alembic.sqlalchemy.org/en/latest/ops.html).
You may also find the other migration scripts under `versions/` useful as a reference.
Finally, Alembic's migration script [tutorial](https://alembic.sqlalchemy.org/en/latest/tutorial.html#create-a-migration-script) may also be useful as a guide.

## Run a migration

From `packages/grz-db`:

```
uv run alembic -c /PATH/TO/alembic.ini upgrade head
```

### Basic alembic.ini

```
[alembic]
# path to migration scripts.
# this is typically a path given in POSIX (e.g. forward slashes)
# format, relative to the token %(here)s which refers to the location of this
# ini file
script_location = .../grz-tools/packages/grz-db/src/grz_db/migrations

# database URL.  This is consumed by the user-maintained env.py script only.
# other means of configuring database URLs may be customized within the env.py
# file.
# See notes in "escaping characters in ini files" for guidelines on
# passwords
sqlalchemy.url = postgresql+psycopg://DBUSER@DBSERVER/DB_NAME
# sqlalchemy.url = sqlite:////PATH/TO/submission.db.sqlite

[loggers]
keys = root,sqlalchemy,alembic

[handlers]
keys = console

[formatters]
keys = generic

[logger_root]
level = WARN
handlers = console
qualname =

[logger_sqlalchemy]
level = WARN
handlers =
qualname = sqlalchemy.engine
propagate = 0

[logger_alembic]
level = INFO
handlers =
qualname = alembic
propagate = 0

[handler_console]
class = StreamHandler
args = (sys.stderr,)
level = NOTSET
formatter = generic

[formatter_generic]
format = %(levelname)-5.5s [%(name)s] %(message)s
```


## General Tips

To easily find the appropriate SQLAlchemy column type for a migration operation, try the following in a REPL after adding the new fields to the particular object:

```py
import sqlmodel.main
from grz_db.models.submission import SubmissionBase
sqlmodel.main.get_column_from_field(SubmissionBase.model_fields["new_column_name"])
```

One can also look at the generated schema for a newly initialized database:

```
sqlite3 submission.db.sqlite .schema
```

For PostgreSQL, Enums must be explicitly created if not added as part of a `create_table` operation (e.g. `add_column`).
See the existing migrations for hints.
