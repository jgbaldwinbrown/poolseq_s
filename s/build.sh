ls *.go | while read i ; do
    go build ${i}
done
