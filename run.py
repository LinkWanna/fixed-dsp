import pathlib
import re

p = pathlib.Path("src/common/tables.rs")
s = p.read_text()
pat = r"pub const SIN_TABLE_I16: \[i16; 513\] = \[(.*?)\n\];"
m = re.search(pat, s, re.S)
assert m, "SIN_TABLE_I16 block not found"
nums = [int(x) for x in re.findall(r"-?\d+", m.group(1))]
assert len(nums) == 513, f"expected 513 items, got {len(nums)}"
lines = ["    " + ", ".join(str(n) for n in nums[i : i + 8]) + "," for i in range(0, len(nums), 8)]
new = "#[rustfmt::skip]\npub const SIN_TABLE_I16: [i16; 513] = [\n" + "\n".join(lines) + "\n];"
s2 = s[: m.start()] + new + s[m.end() :]
p.write_text(s2)
