!/^>/ {
  len=length($0);
  tot = tot+len;
  for (i=len; i>0; --i) {
    ++chars[tolower(substr($0, i, 1))]
  }
}
END {
  print "total:",tot;
  for (c in chars) print c, chars[c];
}
