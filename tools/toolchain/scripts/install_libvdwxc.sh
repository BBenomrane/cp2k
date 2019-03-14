#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"
libvdwxc_ver=${libvdwxc_ver:-master}

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

with_libvdwxc=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_libvdwxc" ] && rm "${BUILDDIR}/setup_libvdwxc"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_libvdwxc" in
    __INSTALL__)
        require_env MPI_CFLAGS
        require_env MPI_LDFLAGS
        require_env MPI_LIBS
        require_env FFTW_LDFLAGS
        require_env FFTW_LIBS
        require_env FFTW_CFLAGS

        echo "==================== Installing libvdwxc ===================="
        pkg_install_dir="${INSTALLDIR}/libvdwxc-${libvdwxc_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if verify_checksums "${install_lock_file}" ; then
            echo "libvdwxc-${libvdwxc_ver} is already installed, skipping it."
        else
            if [ -f libvdwxc-${libvdwxc_ver}-429a80027a2ec2c97e2f6d9a3dc84843f2739865.tar.gz ] ; then
                echo "libvdwxc-${libvdwxc_ver}-429a80027a2ec2c97e2f6d9a3dc84843f2739865.tar.gz is found"
            else
                # do not remove this. They do not publish official version often
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/libvdwxc-master-429a80027a2ec2c97e2f6d9a3dc84843f2739865.tar.gz
            fi

            for patch in "${patches[@]}" ; do
                fname="${patch##*/}"
                if [ -f "${fname}" ] ; then
                    echo "${fname} is found"
                else
                    # parallel build patch
                    download_pkg ${DOWNLOADER_FLAGS} "${patch}"
                fi
            done

            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d libvdwxc-${libvdwxc_ver}-429a80027a2ec2c97e2f6d9a3dc84843f2739865 ] && rm -rf libvdwxc-${libvdwxc_ver}-429a80027a2ec2c97e2f6d9a3dc84843f2739865
            tar -xzf libvdwxc-${libvdwxc_ver}-429a80027a2ec2c97e2f6d9a3dc84843f2739865.tar.gz
            cd libvdwxc-${libvdwxc_ver}-429a80027a2ec2c97e2f6d9a3dc84843f2739865

            for patch in "${patches[@]}" ; do
                patch -p1 < ../"${patch##*/}"
            done
            export
            ./autogen.sh
            ./configure --prefix="${pkg_install_dir}" --with-fftw3
            make -j
            make install
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/$(basename ${SCRIPT_NAME})"
        fi

        LIBVDWXC_CFLAGS="-I${pkg_install_dir}/include"
        LIBVDWXC_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding libvdwxc from system paths ===================="
        check_command pkg-config --modversion libvdwxc
        add_include_from_paths LIBVDWXC_CFLAGS "vdwxc.h" $INCLUDE_PATHS
        add_lib_from_paths LIBVDWXC_LDFLAGS "libvdwxc*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking libvdwxc to user paths ===================="
        pkg_install_dir="$with_libvdwxc"
        check_dir "$pkg_install_dir/lib"
        check_dir "$pkg_install_dir/lib64"
        check_dir "$pkg_install_dir/include"
        LIBVDWXC_CFLAGS="-I'${pkg_install_dir}/include'"
        LIBVDWXC_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_libvdwxc" != "__DONTUSE__" ] ; then
    LIBVDWXC_LIBS="-lvdwxc"
    if [ "$with_libvdwxc" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_libvdwxc"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
    fi
       cat <<EOF > "${BUILDDIR}/setup_libvdwxc"
export LIBVDWXC_CFLAGS="-I$pkg_install_dir/include ${LIBVDWXC_CFLAGS}"
export LIBVDWXC_LDFLAGS="${LIBVDWXC_LDFLAGS}"
export LIBVDWXC_LIBS="${LIBVDWXC_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__LIBVDWXC"
export CP_CFLAGS="\${CP_CFLAGS} ${LIBVDWXC_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBVDWXC_LDFLAGS}"
export CP_LIBS="${LIBVDWXC_LIBS} \${CP_LIBS}"
EOF
        cat "${BUILDDIR}/setup_libvdwxc" >> $SETUPFILE
fi
cd "${ROOTDIR}"
