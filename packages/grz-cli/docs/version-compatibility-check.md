# Version compatibility check

When you run `grz-cli`, it automatically checks whether your installed version is compatible with the current requirements. If your version is too old, too new, or simply behind the recommended release, you will be informed before any submission processing tasks begin.

## What happens when you run grz-cli

On startup, the tool fetches a version policy from a central location and compares it against your installed version. Depending on the result, one of four things will happen:

| Situation | What you'll see | Does execution continue? |
|---|---|---|
| Your version is too old | Error message | No |
| Your version works but is behind the recommendation | Warning message | Yes |
| Your version is newer than the supported maximum | Error message | No |
| Your version is within the supported range | Confirmation message | Yes |

If no policy is currently in effect, the check is silently skipped and execution continues normally.

## What the messages mean

**"Your grz-cli version is not supported"**

Your installed version is below the minimum required. You must upgrade before you can continue.

```
Example message: Your grz-cli version (1.2.0) is not supported. Minimum required version is 1.5.0.
```

**"Upgrading is strongly recommended"**

Your version meets the minimum requirement but is behind the recommended release. You can continue, but upgrading is advised.

```
Example message: You are using grz-cli 1.5.0, but the recommended version is 1.7.0.
Upgrading is strongly recommended.
```

**"grz-cli is within the supported and tested range"**

Your version is fully compatible. No action is needed.

```
Example message: grz-cli 1.7.0 is within the supported and tested range.
```

**"Version is newer than the maximum supported version"**

Your version is ahead of what has been tested and approved. You must downgrade before you can continue.

```
Example message: grz-cli version 2.1.0 is newer than the maximum supported version (2.0.0).
```


## Why is there a maximum version limit?

Most tools only enforce a minimum version. `grz-cli` also enforces an upper bound to prevent untested releases from being used in regulated submission workflows, where the exact toolchain version may need to match what was validated. If you have installed a pre-release or a newer version than your organisation currently supports, downgrade to the approved version before submitting.