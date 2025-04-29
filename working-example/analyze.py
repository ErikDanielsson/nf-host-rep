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
        # Generate character diff
        diff = list(difflib.ndiff(prev, curr))
        
        highlighted_line = ''
        for d in diff:
            if d.startswith('  '):
                highlighted_line += d[-1]
            elif d.startswith('+ '):
                highlighted_line += f"{GREEN}{d[-1]}{RESET}"
            elif d.startswith('- '):
                # You can choose to skip deleted characters, or show them in red
                highlighted_line += f"{RED}{d[-1]}{RESET}"  # Show deletions
                # pass  # Or ignore deletions (just focus on current line)
        
        print(f"Line {i+1}: {highlighted_line}")
        print()
