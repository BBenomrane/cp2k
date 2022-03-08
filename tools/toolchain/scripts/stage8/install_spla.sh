#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

spla_ver="1.5.3"
spla_sha256="527c06e316ce46ec87309a16bfa4138b1abad23fd276fe789c78a2de84f05637"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_SpLA" ] && rm "${BUILDDIR}/setup_SpLA"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_spla" in
  __INSTALL__)
    echo "==================== Installing spla ===================="
    pkg_install_dir="${INSTALLDIR}/SpLA-${spla_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "SpLA-${spla_ver} is already installed, skipping it."
    else
      if [ -f SpLA-${spla_ver}.tar.gz ]; then
        echo "SpLA-${spla_ver}.tar.gz is found"
      else
        download_pkg ${DOWNLOADER_FLAGS} ${spla_sha256} \
          "https://github.com/eth-cscs/Spla/archive/v${spla_ver}.tar.gz" \
          -o SpLA-${spla_ver}.tar.gz

      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d SpLA-${spla_ver} ] && rm -rf SpLA-${spla_ver}
      tar -xzf SpLA-${spla_ver}.tar.gz
      cd spla-${spla_ver}
      mkdir -p build-cpu
      cd build-cpu
      case "${MATH_MODE}" in
        cray)
          EXTRA_CMAKE_FLAGS="-DSPLA_HOST_BLAS=CRAY_LIBSCI"
          ;;
        mkl)
          EXTRA_CMAKE_FLAGS="-DSPLA_HOST_BLAS=MKL"
          ;;
        *)
          EXTRA_CMAKE_FLAGS="-DSPLA_HOST_BLAS=AUTO"
          ;;
      esac
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
        -DSPLA_OMP=OFF \
        -DSPLA_FORTRAN=ON \
        -DSPLA_INSTALL=ON \
        -DSPLA_STATIC=ON \
        ${EXTRA_CMAKE_FLAGS} .. \
        > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make -j $(get_nprocs) install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ..

      if [ "$ENABLE_CUDA" = "__TRUE__" ]; then
        [ -d build-cuda ] && rm -rf "build-cuda"
        mkdir build-cuda
        cd build-cuda
        cmake \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_VERBOSE_MAKEFILE=ON \
          -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
          -DSPLA_OMP=OFF \
          -DSPLA_FORTRAN=ON \
          -DSPLA_INSTALL=ON \
          -DSPLA_STATIC=ON \
          -DSPLA_GPU_BACKEND=CUDA \
          ${EXTRA_CMAKE_FLAGS} .. \
          > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
        make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
        install -d ${pkg_install_dir}/lib/cuda
        [ -f src/libspla.a ] && install -m 644 src/*.a ${pkg_install_dir}/lib/cuda >> install.log 2>&1
        [ -f src/libspla.so ] && install -m 644 src/*.so ${pkg_install_dir}/lib/cuda >> install.log 2>&1
        CP_DFLAGS+=' -D__OFFLOAD_GEMM'
      fi

      if [ "$ENABLE_HIP" = "__TRUE__" ]; then
        case "${GPUVER}" in
          K20X | K40 | K80 | P100 | V100 | A100)
            [ -d build-cuda ] && rm -rf "build-cuda"
            mkdir build-cuda
            cd build-cuda
            cmake \
              -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
              -DCMAKE_INSTALL_LIBDIR=lib \
              -DCMAKE_VERBOSE_MAKEFILE=ON \
              -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
              -DSPLA_OMP=OFF \
              -DSPLA_FORTRAN=ON \
              -DSPLA_INSTALL=ON \
              -DSPLA_STATIC=ON \
              -DSPLA_GPU_BACKEND=ROCM \
              ${EXTRA_CMAKE_FLAGS} .. \
              > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
            make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
            install -d ${pkg_install_dir}/lib/cuda
            [ -f src/libspla.a ] && install -m 644 src/*.a ${pkg_install_dir}/lib/cuda >> install.log 2>&1
            [ -f src/libspla.so ] && install -m 644 src/*.so ${pkg_install_dir}/lib/cuda >> install.log 2>&1
            ;;
          Mi50 | Mi100 | Mi200)
            [ -d build-hip ] && rm -rf "build-hip"
            mkdir build-hip
            cd build-hip
            cmake \
              -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
              -DCMAKE_INSTALL_LIBDIR=lib \
              -DCMAKE_VERBOSE_MAKEFILE=ON \
              -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
              -DSPLA_OMP=OFF \
              -DSPLA_FORTRAN=ON \
              -DSPLA_INSTALL=ON \
              -DSPLA_STATIC=ON \
              -DSPLA_GPU_BACKEND=ROCM \
              ${EXTRA_CMAKE_FLAGS} .. \
              > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
            make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
            install -d ${pkg_install_dir}/lib/rocm
            [ -f src/libspla.a ] && install -m 644 src/*.a ${pkg_install_dir}/lib/rocm >> install.log 2>&1
            [ -f src/libspla.so ] && install -m 644 src/*.so ${pkg_install_dir}/lib/rocm >> install.log 2>&1
            ;;
          *) ;;
        esac
        CP_DFLAGS+=' -D__OFFLOAD_GEMM'
      fi

      # https://github.com/eth-cscs/spla/issues/17
      [ -d "${pkg_install_dir}/lib/cmake/spla/modules" ] && mv "${pkg_install_dir}/lib/cmake/spla/modules" "${pkg_install_dir}/lib/cmake/SPLA/modules"

      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
    fi
    SPLA_ROOT="${pkg_install_dir}"
    SPLA_CFLAGS="-I'${pkg_install_dir}/include/spla'"
    SPLA_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    SPLA_CUDA_LDFLAGS="-L'${pkg_install_dir}/lib/cuda' -Wl,-rpath='${pkg_install_dir}/lib/cuda'"
    SPLA_HIP_LDFLAGS="-L'${pkg_install_dir}/lib/rocm' -Wl,-rpath='${pkg_install_dir}/lib/rocm'"

    ;;
  __SYSTEM__)
    echo "==================== Finding spla from system paths ===================="
    check_command pkg-config --modversion spla
    add_include_from_paths SPLA_CFLAGS "spla.h" $INCLUDE_PATHS
    add_lib_from_paths SPLA_LDFLAGS "libspla.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking spla to user paths ===================="
    pkg_install_dir="$with_spla"

    # use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
    SPLA_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && SPLA_LIBDIR="${pkg_install_dir}/lib64"

    check_dir "${SPLA_LIBDIR}"
    check_dir "${pkg_install_dir}/include/spla"
    SPLA_CFLAGS="-I'${pkg_install_dir}/include/spla'"
    SPLA_LDFLAGS="-L'${SPLA_LIBDIR}' -Wl,-rpath='${SPLA_LIBDIR}'"
    ;;
esac
if [ "$with_spla" != "__DONTUSE__" ]; then
  SPLA_LIBS="-lspla"
  if [ "$with_spla" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_spla"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include/spla"
export SPLA_INCLUDE_DIR="$pkg_install_dir/include/spla"
export SPLA_LIBS="-lspla"
export SPLA_ROOT="${pkg_install_dir}"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}/lib/cmake"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_spla"
export SPLA_CFLAGS="${SPLA_CFLAGS}"
export SPLA_LDFLAGS="${SPLA_LDFLAGS}"
export SPLA_CUDA_LDFLAGS="${SPLA_CUDA_LDFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__SPLA|)"
export CP_CFLAGS="\${CP_CFLAGS} ${SPLA_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_CUDA(${SPLA_CUDA_LDFLAGS}|${SPLA_LDFLAGS})"
export SPLA_LIBRARY="-lspla"
export SPLA_ROOT="$pkg_install_dir"
export SPLA_INCLUDE_DIR="$pkg_install_dir/include/spla"
export SPLA_VERSION=${spla-ver}
export CP_LIBS="IF_MPI(${SPLA_LIBS}|) \${CP_LIBS}"
EOF
  cat "${BUILDDIR}/setup_spla" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_spla"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "spla"
