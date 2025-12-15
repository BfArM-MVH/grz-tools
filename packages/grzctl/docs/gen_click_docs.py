from io import TextIOBase

import click
import grzctl.cli
import mkdocs_gen_files


def recursive_help(command: click.Command, output: TextIOBase, parent: click.Context | None = None, level: int = 2):
    """
    Generates help pages for a click CLI recursively.

    Credit: https://stackoverflow.com/a/58018765
    """
    command_name = "grzctl" if parent is None else command.name
    context = click.Context(command, info_name=command_name, parent=parent)
    header_prefix = "".join(["#"] * level)
    output.write(f"{header_prefix} {command_name}\n")
    output.write("```\n")
    output.write(command.get_help(context))
    output.write("\n```\n\n")
    subcommands = getattr(command, "commands", {}).values()
    for subcommand in subcommands:
        recursive_help(command=subcommand, output=output, parent=context, level=level + 1)


with mkdocs_gen_files.open("cli.md", "w") as cli_docs_file:
    recursive_help(command=grzctl.cli.build_cli(), output=cli_docs_file)


mkdocs_gen_files.set_edit_path("cli.md", "gen_click_docs.py")
