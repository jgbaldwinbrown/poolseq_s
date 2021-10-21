#!/bin/bash
set -e

cp s_syncfile.R ~/mybin/s_syncfile
#cp s_syncfile_win.R ~/mybin/s_syncfile_win
cp s_syncfile_win_typed.R ~/mybin/s_syncfile_win
cp get_partial_s.R ~/mybin/get_partial_s
cp get_partial_s_repl.R ~/mybin/get_partial_s_repl
cp manhat_plot_s_and_p.py ~/mybin/manhat_plot_s_and_p
cp tukey_s.py ~/mybin/tukey_s
cp filter_goods ~/mybin/filter_goods
cp combine_goods ~/mybin/combine_goods
cp manhat_plot_tukey.py ~/mybin/manhat_plot_tukey
go build find_bad_sync_lines.go && cp find_bad_sync_lines ~/mybin/
( cd ne_multi && go build ne_multi.go )
cp ne_multi/ne_multi ~/mybin/ne_multi
