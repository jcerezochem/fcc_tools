#!/bin/bash

git show -s --format=%ci &>/dev/null 
not_git_repo=$?

if (( $not_git_repo )); then

    if [ -e version.f90 ]; then
        echo ""
        echo "*************************************************"
        echo "Not git and version.f90 exists. Reusing it."
        echo "*************************************************"
        echo ""
    else
cat <<EOF > version.f90
    subroutine print_version(unt)
    
        integer,intent(in) :: unt

        write(unt,'(/,X,A)') "GIT VERSION INFO -- Not Available"
        write(unt,'(X,A,/)') " First config date: $(date)"
        
        return
    end subroutine 
EOF

    fi

else

    git_hash=$(git describe --long --dirty --always)
    git_date=$(git show -s --format=%ci)
    user=$USER
    
    if [ -f version.f90 ]; then
        old_git_hash=$(grep "Commit id  :" version.f90)
        old_git_hash=${old_git_hash##*Commit id  : }
        old_git_hash=${old_git_hash%\"}
    else 
        old_git_hash=""
    fi
    
    
    if [ "$old_git_hash" != "$git_hash" ]; then
        echo ""
        echo "*************************************************"
        echo "Git hash changed from previous compilation: "
        echo " Old: $old_git_hash"
        echo " New: $git_hash"
        echo " A new version.f90 file will be generated"
        echo "*************************************************"
        echo ""
cat <<EOF > version.f90
    subroutine print_version(unt)
    
        integer,intent(in) :: unt

        write(unt,'(/,X,A)') "GIT VERSION INFO"
        write(unt,'(X,A)')   " Commit id  : $git_hash"
        write(unt,'(X,A,/)') " Commit date: $git_date"
        
        return
    end subroutine 
EOF
    
    else 
        echo ""
        echo "*************************************************"
        echo "Git hash did not change from previous compilation"
        echo "*************************************************"
        echo ""
    fi

fi