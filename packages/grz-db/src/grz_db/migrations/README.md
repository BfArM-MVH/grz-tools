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

## General Tips

To easily find the appropriate SQLAlchemy column type for a migration operation, try the following in a REPL:

```py
import sqlmodel.main
from grz_db.models.submission import SubmissionBase
sqlmodel.main.get_column_from_field(SubmissionBase.model_fields["new_column_name"])
```

One can also look at the generated schema for a newly initialized database:

```
sqlite3 submission.db.sqlite .schema
```
