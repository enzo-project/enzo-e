import re
_MATCHER = re.compile(r"^(?P<key>.*):(?P<type>.*)=(?P<value>.*)")

def parse_cmake_cache(path):
    out = {}
    with open(path, 'r') as f:
        for line in f:
            if (line.isspace() or (line[:1] == '#') or (line[:2] == '//')):
                # blank line    or  comment          or  docstring
                continue

            match = _MATCHER.match(line.rstrip())
            if match is None:
                raise RuntimeError(
                    "Something went wrong while parsing the following line "
                    f"of '{path}': {line!r}. Lines that aren't blank and "
                    "don't begin with '#' or '//' are expected to follow the "
                    "structure: 'KEY:TYPE=VALUE'")

            # a line consists of 3 parts: KEY:TYPE=VALUE

            # 1. the key value:
            # - it's allowed to be an empty string
            # - when a key contains a colon, it must be enclosed in quotes
            # - it seems the quotes are part of the key name in all other cases
            key = match.group("key")
            if ':' in key:
                assert key[0] == key[-1] == '"'
                key = key[1:-1]

            # 2. get the type:
            # this is used by cmake-gui to help with formatting.
            # Allowed values include:
            #    - BOOL (NOTE: THIS IS HARD TO COERCE!)
            #    - STRING
            #    - FILEPATH
            #    - PATH
            #    - UNINITIALIZED (a variable defined on the command line)
            #    - STATIC (this seems to be a non-configurable value generated)
            #    - INTERNAL
            # we will ignore this for now

            # 3. get the value (it's allowed to be empty)
            value = match.group("value")

            out[key] = value
    return out
