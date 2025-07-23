# Developer guide to migrations

## Create a new migration

From `packages/grz-db`:

```
uv run alembic revision -m "description of migration"
```

To easily find the appropriate SQLAlchemy column type for the migration, try using:

```py
sqlmodel.main.get_column_from_field(SubmissionBase.model_fields["new_column_name"])
```

One can also look at the generated schema for a newly initialized database:

```
sqlite3 submission.db.sqlite .schema
```
