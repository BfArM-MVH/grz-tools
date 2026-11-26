# Deleting submissions on a BfArM request

This page explains how to delete submissions when BfArM asks the GRZ to.
It is written for the data stewards who do this.

- [When you need this](#when-you-need-this)
- [The steps at a glance](#the-steps-at-a-glance)
- [1. Fill in the request](#1-fill-in-the-request)
- [2. Check what would be deleted (dry run)](#2-check-what-would-be-deleted-dry-run)
- [3. Delete](#3-delete)
- [4. Confirm the deletion to BfArM](#4-confirm-the-deletion-to-bfarm)
- [What is deleted, and what is kept](#what-is-deleted-and-what-is-kept)
- [After restoring a backup: `delete check`](#after-restoring-a-backup-delete-check)
- [After three years: `delete purge-records`](#after-three-years-delete-purge-records)
- [Deletions done by hand before `grzctl delete` existed: `delete mark-done`](#deletions-done-by-hand-before-grzctl-delete-existed-delete-mark-done)
- [The quarterly report](#the-quarterly-report)
- [Command reference](#command-reference)

## When you need this

A patient can withdraw their consent, or an LE can find that it sent wrong data.
The LE then fills in BfArM's request form, *Antrag auf Löschung von Daten*.
BfArM sends it to the GRZ as the attachment of an encrypted (S/MIME) e-mail from `PlattformMvGenomseq@bfarm.de`, and tells the GRZ to delete the data.

`grzctl delete` carries out such a request.
It follows BfArM's SOP *Vorgehen zum Löschen von Datensätzen im genomDE-Modellvorhaben während Phase 0* (version 1.0) and the request form in version 1.1.

`grzctl delete` only deletes submissions that are in the submission database.
A submission that is not in the database yet, for example one still waiting in the inbox where the LE noticed an error before processing and vlaidation was carried out, you can delete by hand without involving BfArM.

> `grzctl db submission change-request … Delete` only records that a deletion was requested.
> It deletes nothing.
> Use `grzctl delete` to delete.

## The steps at a glance

| Step | What you do |
| ---- | ----------- |
| 1    | Open BfArM's e-mail in your mail program and save the attached request form (PDF). Copy it to the server where you run grzctl. |
| 2    | Run `grzctl delete template > antrag.yaml` and fill in `antrag.yaml` from the PDF. |
| 3    | Run `grzctl delete run --data-file antrag.yaml --raw-content antrag.pdf --dry-run` and check what would be deleted. |
| 4    | Run the same command without `--dry-run`. This deletes. |
| 5    | Delete the submissions from the GRZ's backups, for example tape backups, as your site's SOP says. grzctl cannot do this. |
| 6    | Run `grzctl delete confirmation --data-file antrag.yaml` and send the text, encrypted, to `PlattformMvGenomseq@bfarm.de`. |
| 7    | Delete `antrag.pdf` and `antrag.yaml` from the server where you ran grzctl. The PDF is kept in the database, and the e-mails stay in the mailbox. |

Your site's SOP says how you delete data from the backups and where you keep the e-mails.

## 1. Fill in the request

grzctl cannot read the PDF.
You copy its values into a YAML file, a plain text file with one value per line.
Get an empty one with:

```console
$ grzctl delete template > antrag.yaml
```

Open `antrag.yaml` in a text editor and replace each `<FILL IN …>` and each `YYYY-MM-DD`:

| Field in `antrag.yaml` | Where you find it |
| ---------------------- | ----------------- |
| `instructed_at`        | The date of BfArM's e-mail. This is the only value that is not on the form. |
| `submitter_id`         | Form: *Meldender LE*. Only the 9-digit number (Institutionskennzeichen), not the place. |
| `requester_name`       | Form: *Kontaktperson*, the name. |
| `requester_email`      | Form: *Kontaktperson*, the e-mail address. |
| `requested_at`         | Form: *Datum des Versands des Antrags auf Löschung*. |
| `reason`               | Form: *Begründung für Löschung*. Copy it as it is. |
| `data_node_id`         | Form: *Beteiligtes GRZ*. Only the ID, for example `GRZM00006`. This is always your own GRZ. |
| `tan_g`                | Form: *tanG des betroffenen genomischen Datensatzes*, or the table the LE attached. One line per tanG. |
| `comment`              | Form: *Kommentare*. You can leave it empty (`""`). |

Leave out the *tanC* and the *Beteiligter KDK* from the form.
They concern the clinical data node, not the GRZ.

A filled-in file looks like this:

```yaml
instructed_at: 2026-11-21
submitter_id: "260914050"
requester_name: "Dr. Erika Mustermann"
requester_email: "erika.mustermann@example.org"
requested_at: 2026-11-20
reason: |
  Widerruf der Einwilligung zur Teilnahme am Modellvorhaben nach § 64e SGB V
data_node_id: "GRZM00006"
tan_g:
  - "3f9a000000000000000000000000000000000000000000000000000000000c21"
  - "8b1e1111111111111111111111111111111111111111111111111111111100d4"
comment: ""
```

## 2. Check what would be deleted (dry run)

A dry run checks everything, but deletes nothing:

```console
$ grzctl delete run --data-file antrag.yaml --raw-content antrag.pdf --dry-run
```

First it checks the file.
If something is missing or wrong, it lists every problem at once, for example a tanG with 63 characters or a date that is not written as `YYYY-MM-DD`.
Fix them and run it again.

Then it looks up the submission of each tanG and shows a table:

```text
Deletion request of 260914050 (LMU München), sent 2026-11-20, instructed 2026-11-21
 tanG       Submission                     State     Result
 3f9a…0c21  260914050_2025-09-15_c64603a7  Reported  will be deleted
 8b1e…00d4  —                              —         not present at GRZM00006
```

The tanGs are shortened on screen, so that they do not end up in terminal logs.

The *Result* column says what will happen:

| Result                                 | Meaning |
| -------------------------------------- | ------- |
| will be deleted                        | The submission is deleted. |
| deletion interrupted, will be finished | An earlier run stopped partway through. This run finishes it. |
| already deleted                        | The submission was deleted before. Nothing happens. |
| not present                            | No submission has this tanG. Check the tanG against the PDF for a typo. If it is right, the GRZ does not have the data, and the confirmation says so. |

A second table shows, for each submission, what each storage place still holds of it, for example `412 objects, 38.2 GB` in the consented archive.

The dry run also warns you:

- **Other submissions of the same case.**
  For example, the request names the first submission of a patient, but not a follow-up.
  grzctl only deletes what the request names.
  If something seems to be missing, ask BfArM.
- **A request that must not be carried out.**
  grzctl stops if the request is addressed to another GRZ, or if a tanG belongs to a submission of another LE.

## 3. Delete

Run the same command without `--dry-run`:

```console
$ grzctl delete run --data-file antrag.yaml --raw-content antrag.pdf
```

grzctl shows the same tables and asks once: *Delete 2 submission(s) for good?*
Answer `y` to delete.
`--yes-i-really-mean-it` skips the question.

For each submission, grzctl:

1. stores the request with the PDF in the database,
2. marks the submission as *Deleting*,
3. deletes its files from the archive and its inbox,
4. checks that nothing is left,
5. removes the personal data from the database,
6. marks the submission as *Deleted*.

**If something goes wrong**, grzctl stops and tells you which submission failed and why.
That submission stays *Deleting*.
Run the same command again: it continues where it stopped, and it does not store the request twice.
Submissions that were deleted already stay deleted.

**Object lock.**
Some storage can lock files so that nobody can delete them.
If a storage place has this switched on, grzctl stops before it deletes anything.
Ask whoever runs your storage.

## 4. Confirm the deletion to BfArM

BfArM needs a confirmation by e-mail that names each tanG.
grzctl writes the text:

```console
$ grzctl delete confirmation --data-file antrag.yaml
```

It looks like this:

```text
Betreff: Löschung von Datensätzen im genomDE MV

Sehr geehrte Damen und Herren,

hiermit bestätigen wir die Löschung der Datensätze gemäß Ihrer Anweisung vom 21.11.2026 (Antrag auf Löschung des Leistungserbringers 260914050 vom 20.11.2026):

  3f9a000000000000000000000000000000000000000000000000000000000c21  gelöscht am 25.11.2026
  8b1e1111111111111111111111111111111111111111111111111111111100d4  liegt beim GRZ GRZM00006 nicht vor

Gelöscht wurden beim GRZ GRZM00006 die genomischen Daten und Metadaten in Archiv und Inbox sowie die personenbezogenen Daten in der Datenbank. Ausgenommen sind die technischen Log-Dateien der Datenübermittlung.

Die Daten wurden ebenso in den Datensicherungen des GRZ GRZM00006 (z. B. Bandsicherungen) gemäß unserer Standort-SOP gelöscht.

Die Löschung wird im Quartalsbericht für Q4/2026 in der Spalte number_of_deletions aufgeführt. Die tanGs werden zur Dokumentation bis längstens 25.11.2029 aufbewahrt.
```

The text says that the data is deleted from the GRZ's backups too.
grzctl cannot do that or check it, so **send the confirmation only after you have deleted the submissions from the backups** (step 5).

Copy the text into an e-mail to `PlattformMvGenomseq@bfarm.de`, add your signature, and **send it encrypted**: it contains the full tanGs.

You can print the text again later, as long as the request is less than three years old.
grzctl refuses to write it while a submission of the request is not deleted yet.

## What is deleted, and what is kept

| Where    | Deleted                                                                                                                                                                                                       | Kept |
| -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---- |
| Archive  | All files of the submission, including `metadata.json`, and older versions of them                                                                                                                            | The log files of the data transfer, which BfArM's SOP allows. An empty file `deleted`, the *tombstone*, which marks that the submission was deleted here. |
| Inbox    | All files of the submission                                                                                                                                                                                   | An empty `metadata.json`, so that the same submission cannot be uploaded again, and the marker `cleaned` |
| Database | The metadata, the local case ID, the donor pseudonyms, the link to the case (the case too, if no other submission belongs to it), the place in the detailed QC queue, the error messages in its state history | The tanG and the request, for three years. The numbers the reports need: type, dates, size, QC results. |

The local case ID becomes `deleted-<submission ID>`, and each donor pseudonym becomes `deleted-<submission ID>-1`, `-2`, and so on.
This way the reports still count each deleted submission, but cannot connect it to other submissions of the same person.

grzctl only reaches the archives, inboxes and database in its configuration.
You delete the submissions from the backups, for example tape backups, and from any other copies by hand, as your site's SOP says, before you send the confirmation.
BfArM's SOP allows at most 90 days for the backups.
Until the data is gone from the backups, nobody may use it, and if a restore brings it back, delete it again with [`delete check`](#after-restoring-a-backup-delete-check).

## After restoring a backup: `delete check`

A backup can bring deleted data back.
BfArM's SOP says you must then delete it again straight away.
So run this after every restore of the database or of an archive or inbox:

```console
$ grzctl delete check
```

It looks for four things:

| Finding              | Meaning                                                                                                                             | `--repair` |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------------------- | ---------- |
| storage restored     | A deleted submission has files in an archive or inbox again.                                                                        | deletes them again |
| database restored    | An archive has a tombstone, but the database does not know the submission was deleted. The database was restored to an older state. | deletes the submission again and removes its personal data again |
| deletion interrupted | A deletion stopped partway through.                                                                                                 | not repaired: run `grzctl delete run` with its request again |
| not in the database  | An archive has a tombstone of a submission the database does not know.                                                              | not repaired: check by hand |

If it finds nothing, it says *Every deleted submission is deleted entirely.*
To delete again what it found:

```console
$ grzctl delete check --repair
```

One gap: a submission that was only in the inbox and never archived leaves no tombstone.
After a database restore, `delete check` cannot find it.

## After three years: `delete purge-records`

grzctl keeps the tanG and the request of each deletion, because BfArM's SOP asks the GRZ to document the deletion.
The same SOP allows this for at most three years.
After that, the records must go.

Run this when you make the quarterly report, which also reminds you of it:

```console
$ grzctl delete purge-records --due-before 2027-06-30 --dry-run
$ grzctl delete purge-records --due-before 2027-06-30
```

`--due-before` takes the end of the quarter after the current one.
That way nothing stays longer than three years until you run it next time.
It deletes a little early rather than late.

grzctl then lists the requests whose e-mails you must delete from the mailbox as well, by the date of BfArM's e-mail and the LE.
grzctl cannot reach the mailbox.

## Deletions done by hand before `grzctl delete` existed: `delete mark-done`

Before `grzctl delete`, deletions were done by hand and recorded with `grzctl db submission change-request … Delete`.
The quarterly report now only counts a deletion that grzctl knows was carried out.
Record each of these old deletions once, with the day it was done:

```console
$ grzctl delete mark-done 260914050_2025-03-01_1a2b3c4d --deleted-on 2026-03-02
```

grzctl refuses if the submission has no *Delete* change request, or if its archive or inbox still holds files.
In that case, delete it with `grzctl delete run` instead.
`mark-done` also removes the submission's personal data from the database.

The date matters: the quarterly report counts the deletion in its quarter, and the three years for the records start on it.

## The quarterly report

`grzctl report quarterly` counts a deletion in `number_of_deletions` of `1-Gesamtübersicht_…tsv`:

- only once the submission is deleted (*Deleted*), not when the deletion is requested;
- in the quarter of the deletion, not in the quarter of the upload.

BfArM's SOP asks for the deletion in the next report after it.

The quarterly report also warns you about:

- submissions with a *Delete* change request that were never deleted with `grzctl delete`;
- records that must be purged before the next report, with the `delete purge-records` command to run.

## Command reference

| Command                                                                                     | What it does |
| ------------------------------------------------------------------------------------------- | ------------ |
| `grzctl delete template`                                                                    | Prints an empty request to fill in. Needs no configuration. |
| `grzctl delete run --data-file FILE --raw-content PDF [--dry-run] [--yes-i-really-mean-it]` | Checks the request and deletes its submissions. `--dry-run` only shows what would be deleted. |
| `grzctl delete confirmation --data-file FILE`                                               | Prints the confirmation e-mail for BfArM. |
| `grzctl delete check [--repair]`                                                            | Finds deleted data that a restore brought back. `--repair` deletes it again. |
| `grzctl delete purge-records --due-before DATE [--dry-run] [--yes-i-really-mean-it]`        | Deletes the records of deletions whose three years end by `DATE`. |
| `grzctl delete mark-done SUBMISSION_ID --deleted-on DATE`                                   | Records a deletion done by hand before `grzctl delete` existed. |

All commands except `template` need the grzctl configuration (`--config`).
