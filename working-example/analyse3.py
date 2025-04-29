import sys
import difflib

GREEN = "\033[92m"
RED = "\033[91m"
RESET = "\033[0m"

if len(sys.argv) != 2:
    print("Usage: python3 highlight_changes.py <filename>")
    sys.exit(1)

filename = sys.argv[1]

with open(filename) as f:
    lines = [line.rstrip('\n') for line in f]

for i in range(1, len(lines)):
    prev, curr = lines[i - 1], lines[i]
    if prev != curr:
        sm = difflib.SequenceMatcher(None, prev, curr)
        old_line = ''
        new_line = ''

        for tag, i1, i2, j1, j2 in sm.get_opcodes():
            if tag == 'equal':
                old_line += prev[i1:i2]
                new_line += curr[j1:j2]
            elif tag == 'delete':
                old_line += f"{RED}{prev[i1:i2]}{RESET}"
            elif tag == 'insert':
                new_line += f"{GREEN}{curr[j1:j2]}{RESET}"
            elif tag == 'replace':
                old_line += f"{RED}{prev[i1:i2]}{RESET}"
                new_line += f"{GREEN}{curr[j1:j2]}{RESET}"

        print(f"\nLine {i}:")
        print("  " + old_line)
        print("  " + new_line)
