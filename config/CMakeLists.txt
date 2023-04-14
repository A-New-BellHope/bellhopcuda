# bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP / BELLHOP3D underwater acoustics simulator
# Copyright (C) 2021-2023 The Regents of the University of California
# Marine Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu
# Based on BELLHOP / BELLHOP3D, which is Copyright (C) 1983-2022 Michael B. Porter
# 
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

# Ensure the clang-format git hook is installed.
set(copy_hook 0)
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/.git/hooks/pre-commit")
    message(WARNING "clang-format pre-commit hook not found, installing")
    set(copy_hook 1)
else()
    file(MD5 "${PROJECT_SOURCE_DIR}/.git/hooks/pre-commit" precommit_md5)
    file(MD5 "${PROJECT_SOURCE_DIR}/clang-format.hook" curhook_md5)
    if(NOT precommit_md5 EQUAL curhook_md5)
        message(WARNING "clang-format pre-commit hook is outdated or a different "
        "pre-commit hook was installed, overwriting with current clang-format hook")
        set(copy_hook 1)
    endif()
endif()
if(copy_hook)
    file(COPY "${PROJECT_SOURCE_DIR}/clang-format.hook"
        DESTINATION "${PROJECT_SOURCE_DIR}/.git/hooks"
    )
    file(RENAME "${PROJECT_SOURCE_DIR}/.git/hooks/clang-format.hook"
        "${PROJECT_SOURCE_DIR}/.git/hooks/pre-commit"
    )
endif()

add_subdirectory(cxx)

if(BHC_ENABLE_CUDA AND NOT(APPLE OR DEFINED ENV{BHC_NO_CUDA}))
    add_subdirectory(cuda)
endif()

if(BHC_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
