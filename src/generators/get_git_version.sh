#!/bin/bash

# git_hash=$(git rev-parse HEAD)
git_hash=$(git describe --long --dirty --always)
git_date=$(git show -s --format=%ci)
user=$USER

#echo "! Automatically generated subroutine with git hash" > version.f90
#echo "" >> version.f90

#echo "module version" >> version.f90

#echo "    contains" >> version.f90
#echo "    " >> version.f90
cat <<EOF > version.f90
    subroutine print_version()

        write(0,'(/A)') "GIT VERSION INFO"
        write(0,'(A)')   " Commit id  : $git_hash"
        write(0,'(A)')   " Commit date: $git_date"
        
        return
    end subroutine 
EOF

#        write(6,'(13X,A)')   "Compilation date: $date"
#        write(6,'(13X,A)')   "Flags: $@"

#echo "end module" >> version.f90
