---
name: flex-testers-stay-local
description: the user does not want the flex unit testers (e.g. simple_flex_cls_split_tester.f90 and its suite registration) pushed to origin; they stay on the local branch
metadata:
  type: feedback
---

The flex testers (`simple_flex_cls_split_tester.f90`, the `add_suite` / `use` registration lines in
`simple_commanders_test_class.f90` and the strings in `simple_test_ui_class.f90`) are NOT part of what
gets pushed to origin. User, 2026-10-05: "did I not tell you to not push the tester?"

**Why:** the suites are the user's development instrument on the branch; master gets the method, the
strategy, the commander/UI entry and the notes.

**How to apply:** when squashing a flex branch for origin, exclude the tester file and the two
registration-line files from the commit; keep them on the local branch. Same rule as
[[never-push-to-origin-without-permission]]: list the files before any push.
