#!/usr/bin/env python3
import sys

# C++を圧縮

RAW_PREFIXES = ["u8R\"", "uR\"", "UR\"", "LR\"", "R\""]
STR_PREFIXES = ["u8\"", "u\"", "U\"", "L\"", "\""]
CHAR_PREFIXES = ["u'", "U'", "L'", "'"]

PUNCTUATORS = [
    "%:%:", "<<=", ">>=", "<=>", "->*", "...",
    "::", "->", "++", "--", "<<", ">>", "<=", ">=", "==", "!=",
    "&&", "||", "*=", "/=", "%=", "+=", "-=", "&=", "^=", "|=",
    ".*", "<:", ":>", "<%", "%>", "%:", "##",
]

MERGE_PUNCT = set(PUNCTUATORS)
MERGE_PUNCT.add("//")
MERGE_PUNCT.add("/*")


def parse_raw_string(s: str, i: int, prefix_len: int):
    j = i + prefix_len
    k = s.find("(", j)
    if k == -1:
        return None

    delim = s[j:k]
    end_token = ")" + delim + "\""
    end = s.find(end_token, k + 1)
    if end == -1:
        return len(s)

    return end + len(end_token)


def parse_escaped_literal(s: str, i: int, prefix: str):
    j = i + len(prefix)
    quote = prefix[-1]

    while j < len(s):
        if s[j] == "\\":
            j += 2
        elif s[j] == quote:
            return j + 1
        else:
            j += 1

    return len(s)


def parse_literal(s: str, i: int):
    for p in RAW_PREFIXES:
        if s.startswith(p, i):
            return ("string", parse_raw_string(s, i, len(p)))

    for p in STR_PREFIXES:
        if s.startswith(p, i):
            return ("string", parse_escaped_literal(s, i, p))

    for p in CHAR_PREFIXES:
        if s.startswith(p, i):
            return ("char", parse_escaped_literal(s, i, p))

    return None


def remove_comments(src: str) -> str:
    src = src.replace("\r\n", "\n").replace("\r", "\n")
    out = []
    i = 0
    n = len(src)

    while i < n:
        lit = parse_literal(src, i)
        if lit is not None:
            _, end = lit
            out.append(src[i:end])
            i = end
            continue

        if i + 1 < n and src[i] == "/" and src[i + 1] == "/":
            i += 2
            while i < n and src[i] != "\n":
                i += 1
            continue

        if i + 1 < n and src[i] == "/" and src[i + 1] == "*":
            i += 2
            while i < n:
                if i + 1 < n and src[i] == "*" and src[i + 1] == "/":
                    i += 2
                    break
                if src[i] == "\n":
                    out.append("\n")
                i += 1
            continue

        out.append(src[i])
        i += 1

    return "".join(out)


def is_ident_start(c: str) -> bool:
    return c == "_" or c.isalpha() or ord(c) >= 128


def is_ident_continue(c: str) -> bool:
    return c == "_" or c.isalnum() or ord(c) >= 128


def read_identifier(s: str, i: int) -> int:
    j = i + 1
    while j < len(s) and is_ident_continue(s[j]):
        j += 1
    return j


def read_number(s: str, i: int) -> int:
    j = i
    if s[j] == ".":
        j += 1

    while j < len(s):
        c = s[j]
        if c.isalnum() or c in "._'":
            j += 1
        elif c in "+-" and j > i and s[j - 1] in "eEpP":
            j += 1
        else:
            break

    return j


def read_punctuator(s: str, i: int) -> int:
    for p in PUNCTUATORS:
        if s.startswith(p, i):
            return i + len(p)
    return i + 1


def tokenize_code(s: str):
    i = 0
    n = len(s)

    while i < n:
        if s[i].isspace():
            i += 1
            continue

        lit = parse_literal(s, i)
        if lit is not None:
            kind, end = lit
            yield (kind, s[i:end])
            i = end
            continue

        c = s[i]

        if is_ident_start(c):
            end = read_identifier(s, i)
            yield ("id", s[i:end])
            i = end
            continue

        if c.isdigit() or (c == "." and i + 1 < n and s[i + 1].isdigit()):
            end = read_number(s, i)
            yield ("number", s[i:end])
            i = end
            continue

        end = read_punctuator(s, i)
        yield ("punct", s[i:end])
        i = end


def needs_space(prev, cur) -> bool:
    pk, pt = prev
    ck, ct = cur

    if pk in {"id", "number"} and ck in {"id", "number"}:
        return True

    if pk in {"string", "char", "number"} and ck == "id":
        return True

    if pk == "number" and ct.startswith("."):
        return True

    if pk == "punct" and ck == "punct" and (pt + ct) in MERGE_PUNCT:
        return True

    return False


def minify_code(s: str) -> str:
    out = []
    prev = None

    for tok in tokenize_code(s):
        if prev is not None and needs_space(prev, tok):
            out.append(" ")
        out.append(tok[1])
        prev = tok

    return "".join(out)


def is_pp_start(line: str) -> bool:
    return line.lstrip(" \t\v\f").startswith("#")


def minify_cpp(src: str) -> str:
    src = remove_comments(src)
    lines = src.splitlines(keepends=True)

    out = []
    code_buf = []
    i = 0

    def flush_code():
        nonlocal code_buf
        if not code_buf:
            return
        chunk = "".join(code_buf)
        code_buf = []
        m = minify_code(chunk)
        if m:
            out.append(m)

    while i < len(lines):
        line = lines[i]

        if is_pp_start(line):
            flush_code()
            if out and not out[-1].endswith("\n"):
                out.append("\n")

            first = True
            while True:
                current = lines[i]
                if first:
                    current = current.lstrip(" \t\v\f")
                    first = False
                out.append(current)

                continued = current.endswith("\\\n") or current.endswith("\\")
                i += 1
                if not continued or i >= len(lines):
                    break
            continue

        code_buf.append(line)
        i += 1

    flush_code()

    result = "".join(out)
    if not result.endswith("\n"):
        result += "\n"
    return result


def main():
    src = sys.stdin.read()
    sys.stdout.write(minify_cpp(src))


if __name__ == "__main__":
    main()
