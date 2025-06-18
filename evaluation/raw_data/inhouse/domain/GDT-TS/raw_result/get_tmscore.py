import sys

if len(sys.argv) != 2:
    print(f"Usage: python {sys.argv[0]} <input_file>")
    sys.exit(1)

input_file = sys.argv[1]

with open(input_file, "r") as f:
    for line in f:
        parts = line.strip().split()
        if not parts:
            continue
        target = parts[0]
        try:
            scores = list(map(float, parts[1:]))
            top1 = scores[0]
            best = max(scores)
            print(f"{target}\t{top1:.4f}\t{best:.4f}")
        except ValueError:
            print(f"Skipping line due to parse error: {line.strip()}")

