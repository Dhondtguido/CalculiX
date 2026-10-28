add_library(flags INTERFACE)

target_compile_options(flags INTERFACE

  # ====================
  # GNU (gcc + gfortran)
  # ====================

  $<$<AND:$<BOOL:${CALCULIX_INT64}>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-fdefault-integer-8>

  # Release
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:C,GNU>>:>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-cpp>

  # RelWithDebInfo
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:C,GNU>>:-O3 -fno-omit-frame-pointer>
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-cpp -O3 -fno-omit-frame-pointer>

  # Debug
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:C,GNU>>:-Wall>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-cpp -Wall>


  # ===============================
  # Clang toolchain (clang + flang)
  # ===============================

  $<$<AND:$<BOOL:${CALCULIX_INT64}>,$<COMPILE_LANG_AND_ID:Fortran,LLVMFlang>>:-i8>

  # Release
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:C,Clang,AppleClang>>:>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:Fortran,LLVMFlang>>:-cpp>

  # RelWithDebInfo
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:C,Clang,AppleClang>>:-O3 -fno-omit-frame-pointer>
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:Fortran,LLVMFlang>>:-cpp -O3 -fno-omit-frame-pointer>

  # Debug
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:C,Clang,AppleClang>>:-Wall>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:Fortran,LLVMFlang>>:-cpp -Wall>


  # ======================
  # Intel LLVM (icx + ifx)
  # ======================

  $<$<AND:$<BOOL:${CALCULIX_INT64}>,$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>>:-i8>

  # Release
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:C,IntelLLVM>>:-fpic>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>>:-fpp -fpic>

  # RelWithDebInfo
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:C,IntelLLVM>>:-O3 -fpic -fno-omit-frame-pointer>
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>>:-fpp -O3 -fpic -fno-omit-frame-pointer>

  # Debug
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:C,IntelLLVM>>:-Wall -fpic -fno-omit-frame-pointer -fsanitize=address>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>>:-fpp -Wall -fpic -fno-omit-frame-pointer -fsanitize=address>


  # ===========================
  # Intel classic (icc + ifort)
  # ===========================

  $<$<AND:$<BOOL:${CALCULIX_INT64}>,$<COMPILE_LANG_AND_ID:Fortran,Intel>>:-i8>

  # Release
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:C,Intel>>:>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:Fortran,Intel>>:-fpp>

  # RelWithDebInfo
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:C,Intel>>:-O3 -fno-omit-frame-pointer>
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:Fortran,Intel>>:-fpp -O3 -fno-omit-frame-pointer>

  # Debug
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:C,Intel>>:-Wall>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:Fortran,Intel>>:-fpp -Wall>

)


target_link_options(flags INTERFACE

  # ====================
  # GNU (gcc + gfortran)
  # ====================

  # Release
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:C,GNU>>:>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:>

  # RelWithDebInfo
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:C,GNU>>:>
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:>

  # Debug
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:C,GNU>>:>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:>


  # ===============================
  # Clang toolchain (clang + flang)
  # ===============================

  # Release
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:C,Clang,AppleClang>>:>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:Fortran,LLVMFlang>>:>

  # RelWithDebInfo
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:C,Clang,AppleClang>>:>
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:Fortran,LLVMFlang>>:>

  # Debug
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:C,Clang,AppleClang>>:>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:Fortran,LLVMFlang>>:>


  # ======================
  # Intel LLVM (icx + ifx)
  # ======================

  # Release
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:C,IntelLLVM>>:>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>>:>

  # RelWithDebInfo
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:C,IntelLLVM>>:>
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>>:>

  # Debug
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:C,IntelLLVM>>:-fsanitize=address>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>>:-fsanitize=address>


  # ===========================
  # Intel classic (icc + ifort)
  # ===========================

  # Release
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:C,Intel>>:>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:Fortran,Intel>>:>

  # RelWithDebInfo
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:C,Intel>>:>
  $<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANG_AND_ID:Fortran,Intel>>:>

  # Debug
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:C,Intel>>:>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:Fortran,Intel>>:>

)
