if [ $version ]
then
    export version=${version//[[:alpha:]]/}
    echo "get env version="$version
else
    export version="0.0.0";
    echo "set env version="$version
fi