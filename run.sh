# Script for running Jasmine
if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
    TOTAL_MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    MAX_HEAP="$((TOTAL_MEM_KB * 8 / 10 / 1024))m"
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
    TOTAL_MEM_BYTES=$(sysctl -n hw.memsize)
    MAX_HEAP="$((TOTAL_MEM_BYTES * 8 / 10 / 1024 / 1024))m"
fi

java -Xmx${MAX_HEAP} -cp $BINDIR/src:$BINDIR/Iris/src Main "${@:1}"
