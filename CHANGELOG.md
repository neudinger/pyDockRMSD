# Changelog

All notable changes to this project are documented in this file. On the [releases page](https://github.com/neudinger/pyDockRMSD/releases/) you can see all released versions and download the [latest version](https://github.com/neudinger/pyDockRMSD/releases/latest).

## [1.0.0] - 2022-03-06

**This release change the memory allocation types.**

To process large file, stack allocations must be changed to heap allocation. But memory management in c must be made manually with free.

- Replace all stack allocation to heap allocations.

    Change all allocation with [alloca](https://man7.org/linux/man-pages/man3/alloca.3.html) to [malloc](https://en.cppreference.com/w/c/memory/malloc).

    For more information visit the issue <https://github.com/neudinger/pyDockRMSD/issues/3>.

- Add main in comment for test and debug developments.

- Remove all _alloca

    Delete [_alloca](https://docs.microsoft.com/en-us/cpp/c-runtime-library/reference/alloca?view=msvc-170) to remove compiler dependency. 
