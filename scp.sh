#!/usr/bin/expect

set directory [lindex $argv 0]

spawn scp -r ./$directory rl705350@abel.polytech.unice.fr:~/
expect "password"
send "SUk.w@xe\r"
interact