# SFREE Macro Bug Fix

## Problem
The `SFREE` macro in CalculiX doesn't set the pointer to NULL after freeing memory, which can lead to double-free crashes if `SFREE` is called multiple times on the same variable.

## Fix
Changed line 64 in `src/CalculiX.h`:

**Before:**
```c
#define SFREE(a) u_free(a,__FILE__,__LINE__,#a)
```

**After:**
```c
#define SFREE(a) ((a)=u_free(a,__FILE__,__LINE__,#a))
```

## Files Changed
- `src/CalculiX.h` - SFREE macro fix
- `test_sfree.c` - Test to verify the fix works correctly

## Test
Run `./test_sfree` to verify the fix. It demonstrates:
1. The OLD macro leaks dangling pointers
2. The FIXED macro properly sets pointers to NULL

## Branch
`fix/sfree-macro-bug` in `~/dev/CalculiX`

## To Push and Create PR
```bash
cd ~/dev/CalculiX
git push -u origin fix/sfree-macro-bug
```

Then create a PR from the `fix/sfree-macro-bug` branch to `master`.
