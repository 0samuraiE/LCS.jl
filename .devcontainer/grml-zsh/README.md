# grml-zsh Devcontainer Feature

Configures zsh with the [grml zsh configuration](https://grml.org/zsh/) and popular plugins.

## Features

- **grml base configuration**: Powerful zsh configuration with sensible defaults
- **zsh-completions**: Additional completion definitions
- **zsh-syntax-highlighting**: Fish-like syntax highlighting
- **zsh-autosuggestions**: Fish-like autosuggestions
- **zsh-history-substring-search**: Improved history search

## Usage

```json
{
    "features": {
        "./grml-zsh": {}
    }
}
```

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| enableCompletions | boolean | true | Enable zsh-completions plugin |
| enableSyntaxHighlighting | boolean | true | Enable syntax highlighting |
| enableAutosuggestions | boolean | true | Enable autosuggestions |
| enableHistorySubstring | boolean | true | Enable history substring search |

## Example

Disable specific plugins:

```json
{
    "features": {
        "./grml-zsh": {
            "enableHistorySubstring": false
        }
    }
}
```

## Prerequisites

- zsh must be installed in the base image
- git, wget, and ca-certificates must be available

## Configuration Files

The feature creates three configuration files in the user's home directory:

- `~/.zshrc`: Main grml configuration (downloaded from grml.org)
- `~/.zshrc.pre`: Executed before grml config (PATH, fpath)
- `~/.zshrc.local`: Executed after grml config (plugin loading)

## References

- [grml zsh](https://grml.org/zsh/)
- [zsh-users organization](https://github.com/zsh-users)
