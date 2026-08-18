import { readFileSync } from "node:fs";
import { parentPort } from "node:worker_threads";
var Iovec = class Iovec {
	static read_bytes(view, ptr) {
		const iovec = new Iovec();
		iovec.buf = view.getUint32(ptr, true);
		iovec.buf_len = view.getUint32(ptr + 4, true);
		return iovec;
	}
	static read_bytes_array(view, ptr, len) {
		const iovecs = [];
		for (let i = 0; i < len; i++) iovecs.push(Iovec.read_bytes(view, ptr + 8 * i));
		return iovecs;
	}
};
var Ciovec = class Ciovec {
	static read_bytes(view, ptr) {
		const iovec = new Ciovec();
		iovec.buf = view.getUint32(ptr, true);
		iovec.buf_len = view.getUint32(ptr + 4, true);
		return iovec;
	}
	static read_bytes_array(view, ptr, len) {
		const iovecs = [];
		for (let i = 0; i < len; i++) iovecs.push(Ciovec.read_bytes(view, ptr + 8 * i));
		return iovecs;
	}
};
var Subscription = class Subscription {
	static read_bytes(view, ptr) {
		return new Subscription(view.getBigUint64(ptr, true), view.getUint8(ptr + 8), view.getUint32(ptr + 16, true), view.getBigUint64(ptr + 24, true), view.getUint16(ptr + 36, true));
	}
	constructor(userdata, eventtype, clockid, timeout, flags) {
		this.userdata = userdata;
		this.eventtype = eventtype;
		this.clockid = clockid;
		this.timeout = timeout;
		this.flags = flags;
	}
};
var Event = class {
	write_bytes(view, ptr) {
		view.setBigUint64(ptr, this.userdata, true);
		view.setUint16(ptr + 8, this.error, true);
		view.setUint8(ptr + 10, this.eventtype);
	}
	constructor(userdata, error, eventtype) {
		this.userdata = userdata;
		this.error = error;
		this.eventtype = eventtype;
	}
};
//#endregion
//#region ../swig-flow-refactor/node_modules/@bjorn3/browser_wasi_shim/dist/debug.js
let Debug = class Debug {
	enable(enabled) {
		this.log = createLogger(enabled === void 0 ? true : enabled, this.prefix);
	}
	get enabled() {
		return this.isEnabled;
	}
	constructor(isEnabled) {
		this.isEnabled = isEnabled;
		this.prefix = "wasi:";
		this.enable(isEnabled);
	}
};
function createLogger(enabled, prefix) {
	if (enabled) return console.log.bind(console, "%c%s", "color: #265BA0", prefix);
	else return () => {};
}
const debug = new Debug(false);
//#endregion
//#region ../swig-flow-refactor/node_modules/@bjorn3/browser_wasi_shim/dist/wasi.js
var WASIProcExit = class extends Error {
	constructor(code) {
		super("exit with exit code " + code);
		this.code = code;
	}
};
let WASI = class WASI {
	start(instance) {
		this.inst = instance;
		try {
			instance.exports._start();
			return 0;
		} catch (e) {
			if (e instanceof WASIProcExit) return e.code;
			else throw e;
		}
	}
	initialize(instance) {
		this.inst = instance;
		if (instance.exports._initialize) instance.exports._initialize();
	}
	constructor(args, env, fds, options = {}) {
		this.args = [];
		this.env = [];
		this.fds = [];
		debug.enable(options.debug);
		this.args = args;
		this.env = env;
		this.fds = fds;
		const self = this;
		this.wasiImport = {
			args_sizes_get(argc, argv_buf_size) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				buffer.setUint32(argc, self.args.length, true);
				let buf_size = 0;
				for (const arg of self.args) buf_size += arg.length + 1;
				buffer.setUint32(argv_buf_size, buf_size, true);
				debug.log(buffer.getUint32(argc, true), buffer.getUint32(argv_buf_size, true));
				return 0;
			},
			args_get(argv, argv_buf) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				const orig_argv_buf = argv_buf;
				for (let i = 0; i < self.args.length; i++) {
					buffer.setUint32(argv, argv_buf, true);
					argv += 4;
					const arg = new TextEncoder().encode(self.args[i]);
					buffer8.set(arg, argv_buf);
					buffer.setUint8(argv_buf + arg.length, 0);
					argv_buf += arg.length + 1;
				}
				if (debug.enabled) debug.log(new TextDecoder("utf-8").decode(buffer8.slice(orig_argv_buf, argv_buf)));
				return 0;
			},
			environ_sizes_get(environ_count, environ_size) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				buffer.setUint32(environ_count, self.env.length, true);
				let buf_size = 0;
				for (const environ of self.env) buf_size += new TextEncoder().encode(environ).length + 1;
				buffer.setUint32(environ_size, buf_size, true);
				debug.log(buffer.getUint32(environ_count, true), buffer.getUint32(environ_size, true));
				return 0;
			},
			environ_get(environ, environ_buf) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				const orig_environ_buf = environ_buf;
				for (let i = 0; i < self.env.length; i++) {
					buffer.setUint32(environ, environ_buf, true);
					environ += 4;
					const e = new TextEncoder().encode(self.env[i]);
					buffer8.set(e, environ_buf);
					buffer.setUint8(environ_buf + e.length, 0);
					environ_buf += e.length + 1;
				}
				if (debug.enabled) debug.log(new TextDecoder("utf-8").decode(buffer8.slice(orig_environ_buf, environ_buf)));
				return 0;
			},
			clock_res_get(id, res_ptr) {
				let resolutionValue;
				switch (id) {
					case 1:
						resolutionValue = 5000n;
						break;
					case 0:
						resolutionValue = 1000000n;
						break;
					default: return 52;
				}
				new DataView(self.inst.exports.memory.buffer).setBigUint64(res_ptr, resolutionValue, true);
				return 0;
			},
			clock_time_get(id, precision, time) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				if (id === 0) buffer.setBigUint64(time, BigInt((/* @__PURE__ */ new Date()).getTime()) * 1000000n, true);
				else if (id == 1) {
					let monotonic_time;
					try {
						monotonic_time = BigInt(Math.round(performance.now() * 1e6));
					} catch (e) {
						monotonic_time = 0n;
					}
					buffer.setBigUint64(time, monotonic_time, true);
				} else buffer.setBigUint64(time, 0n, true);
				return 0;
			},
			fd_advise(fd, offset, len, advice) {
				if (self.fds[fd] != void 0) return 0;
				else return 8;
			},
			fd_allocate(fd, offset, len) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_allocate(offset, len);
				else return 8;
			},
			fd_close(fd) {
				if (self.fds[fd] != void 0) {
					const ret = self.fds[fd].fd_close();
					self.fds[fd] = void 0;
					return ret;
				} else return 8;
			},
			fd_datasync(fd) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_sync();
				else return 8;
			},
			fd_fdstat_get(fd, fdstat_ptr) {
				if (self.fds[fd] != void 0) {
					const { ret, fdstat } = self.fds[fd].fd_fdstat_get();
					if (fdstat != null) fdstat.write_bytes(new DataView(self.inst.exports.memory.buffer), fdstat_ptr);
					return ret;
				} else return 8;
			},
			fd_fdstat_set_flags(fd, flags) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_fdstat_set_flags(flags);
				else return 8;
			},
			fd_fdstat_set_rights(fd, fs_rights_base, fs_rights_inheriting) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_fdstat_set_rights(fs_rights_base, fs_rights_inheriting);
				else return 8;
			},
			fd_filestat_get(fd, filestat_ptr) {
				if (self.fds[fd] != void 0) {
					const { ret, filestat } = self.fds[fd].fd_filestat_get();
					if (filestat != null) filestat.write_bytes(new DataView(self.inst.exports.memory.buffer), filestat_ptr);
					return ret;
				} else return 8;
			},
			fd_filestat_set_size(fd, size) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_filestat_set_size(size);
				else return 8;
			},
			fd_filestat_set_times(fd, atim, mtim, fst_flags) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_filestat_set_times(atim, mtim, fst_flags);
				else return 8;
			},
			fd_pread(fd, iovs_ptr, iovs_len, offset, nread_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const iovecs = Iovec.read_bytes_array(buffer, iovs_ptr, iovs_len);
					let nread = 0;
					for (const iovec of iovecs) {
						const { ret, data } = self.fds[fd].fd_pread(iovec.buf_len, offset);
						if (ret != 0) {
							buffer.setUint32(nread_ptr, nread, true);
							return ret;
						}
						buffer8.set(data, iovec.buf);
						nread += data.length;
						offset += BigInt(data.length);
						if (data.length != iovec.buf_len) break;
					}
					buffer.setUint32(nread_ptr, nread, true);
					return 0;
				} else return 8;
			},
			fd_prestat_get(fd, buf_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const { ret, prestat } = self.fds[fd].fd_prestat_get();
					if (prestat != null) prestat.write_bytes(buffer, buf_ptr);
					return ret;
				} else return 8;
			},
			fd_prestat_dir_name(fd, path_ptr, path_len) {
				if (self.fds[fd] != void 0) {
					const { ret, prestat } = self.fds[fd].fd_prestat_get();
					if (prestat == null) return ret;
					const prestat_dir_name = prestat.inner.pr_name;
					new Uint8Array(self.inst.exports.memory.buffer).set(prestat_dir_name.slice(0, path_len), path_ptr);
					return prestat_dir_name.byteLength > path_len ? 37 : 0;
				} else return 8;
			},
			fd_pwrite(fd, iovs_ptr, iovs_len, offset, nwritten_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const iovecs = Ciovec.read_bytes_array(buffer, iovs_ptr, iovs_len);
					let nwritten = 0;
					for (const iovec of iovecs) {
						const data = buffer8.slice(iovec.buf, iovec.buf + iovec.buf_len);
						const { ret, nwritten: nwritten_part } = self.fds[fd].fd_pwrite(data, offset);
						if (ret != 0) {
							buffer.setUint32(nwritten_ptr, nwritten, true);
							return ret;
						}
						nwritten += nwritten_part;
						offset += BigInt(nwritten_part);
						if (nwritten_part != data.byteLength) break;
					}
					buffer.setUint32(nwritten_ptr, nwritten, true);
					return 0;
				} else return 8;
			},
			fd_read(fd, iovs_ptr, iovs_len, nread_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const iovecs = Iovec.read_bytes_array(buffer, iovs_ptr, iovs_len);
					let nread = 0;
					for (const iovec of iovecs) {
						const { ret, data } = self.fds[fd].fd_read(iovec.buf_len);
						if (ret != 0) {
							buffer.setUint32(nread_ptr, nread, true);
							return ret;
						}
						buffer8.set(data, iovec.buf);
						nread += data.length;
						if (data.length != iovec.buf_len) break;
					}
					buffer.setUint32(nread_ptr, nread, true);
					return 0;
				} else return 8;
			},
			fd_readdir(fd, buf, buf_len, cookie, bufused_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					let bufused = 0;
					while (true) {
						const { ret, dirent } = self.fds[fd].fd_readdir_single(cookie);
						if (ret != 0) {
							buffer.setUint32(bufused_ptr, bufused, true);
							return ret;
						}
						if (dirent == null) break;
						if (buf_len - bufused < dirent.head_length()) {
							bufused = buf_len;
							break;
						}
						const head_bytes = new ArrayBuffer(dirent.head_length());
						dirent.write_head_bytes(new DataView(head_bytes), 0);
						buffer8.set(new Uint8Array(head_bytes).slice(0, Math.min(head_bytes.byteLength, buf_len - bufused)), buf);
						buf += dirent.head_length();
						bufused += dirent.head_length();
						if (buf_len - bufused < dirent.name_length()) {
							bufused = buf_len;
							break;
						}
						dirent.write_name_bytes(buffer8, buf, buf_len - bufused);
						buf += dirent.name_length();
						bufused += dirent.name_length();
						cookie = dirent.d_next;
					}
					buffer.setUint32(bufused_ptr, bufused, true);
					return 0;
				} else return 8;
			},
			fd_renumber(fd, to) {
				if (self.fds[fd] != void 0 && self.fds[to] != void 0) {
					const ret = self.fds[to].fd_close();
					if (ret != 0) return ret;
					self.fds[to] = self.fds[fd];
					self.fds[fd] = void 0;
					return 0;
				} else return 8;
			},
			fd_seek(fd, offset, whence, offset_out_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const { ret, offset: offset_out } = self.fds[fd].fd_seek(offset, whence);
					buffer.setBigInt64(offset_out_ptr, offset_out, true);
					return ret;
				} else return 8;
			},
			fd_sync(fd) {
				if (self.fds[fd] != void 0) return self.fds[fd].fd_sync();
				else return 8;
			},
			fd_tell(fd, offset_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const { ret, offset } = self.fds[fd].fd_tell();
					buffer.setBigUint64(offset_ptr, offset, true);
					return ret;
				} else return 8;
			},
			fd_write(fd, iovs_ptr, iovs_len, nwritten_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const iovecs = Ciovec.read_bytes_array(buffer, iovs_ptr, iovs_len);
					let nwritten = 0;
					for (const iovec of iovecs) {
						const data = buffer8.slice(iovec.buf, iovec.buf + iovec.buf_len);
						const { ret, nwritten: nwritten_part } = self.fds[fd].fd_write(data);
						if (ret != 0) {
							buffer.setUint32(nwritten_ptr, nwritten, true);
							return ret;
						}
						nwritten += nwritten_part;
						if (nwritten_part != data.byteLength) break;
					}
					buffer.setUint32(nwritten_ptr, nwritten, true);
					return 0;
				} else return 8;
			},
			path_create_directory(fd, path_ptr, path_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					return self.fds[fd].path_create_directory(path);
				} else return 8;
			},
			path_filestat_get(fd, flags, path_ptr, path_len, filestat_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					const { ret, filestat } = self.fds[fd].path_filestat_get(flags, path);
					if (filestat != null) filestat.write_bytes(buffer, filestat_ptr);
					return ret;
				} else return 8;
			},
			path_filestat_set_times(fd, flags, path_ptr, path_len, atim, mtim, fst_flags) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					return self.fds[fd].path_filestat_set_times(flags, path, atim, mtim, fst_flags);
				} else return 8;
			},
			path_link(old_fd, old_flags, old_path_ptr, old_path_len, new_fd, new_path_ptr, new_path_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[old_fd] != void 0 && self.fds[new_fd] != void 0) {
					const old_path = new TextDecoder("utf-8").decode(buffer8.slice(old_path_ptr, old_path_ptr + old_path_len));
					const new_path = new TextDecoder("utf-8").decode(buffer8.slice(new_path_ptr, new_path_ptr + new_path_len));
					const { ret, inode_obj } = self.fds[old_fd].path_lookup(old_path, old_flags);
					if (inode_obj == null) return ret;
					return self.fds[new_fd].path_link(new_path, inode_obj, false);
				} else return 8;
			},
			path_open(fd, dirflags, path_ptr, path_len, oflags, fs_rights_base, fs_rights_inheriting, fd_flags, opened_fd_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					debug.log(path);
					const { ret, fd_obj } = self.fds[fd].path_open(dirflags, path, oflags, fs_rights_base, fs_rights_inheriting, fd_flags);
					if (ret != 0) return ret;
					self.fds.push(fd_obj);
					const opened_fd = self.fds.length - 1;
					buffer.setUint32(opened_fd_ptr, opened_fd, true);
					return 0;
				} else return 8;
			},
			path_readlink(fd, path_ptr, path_len, buf_ptr, buf_len, nread_ptr) {
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					debug.log(path);
					const { ret, data } = self.fds[fd].path_readlink(path);
					if (data != null) {
						const data_buf = new TextEncoder().encode(data);
						if (data_buf.length > buf_len) {
							buffer.setUint32(nread_ptr, 0, true);
							return 8;
						}
						buffer8.set(data_buf, buf_ptr);
						buffer.setUint32(nread_ptr, data_buf.length, true);
					}
					return ret;
				} else return 8;
			},
			path_remove_directory(fd, path_ptr, path_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					return self.fds[fd].path_remove_directory(path);
				} else return 8;
			},
			path_rename(fd, old_path_ptr, old_path_len, new_fd, new_path_ptr, new_path_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0 && self.fds[new_fd] != void 0) {
					const old_path = new TextDecoder("utf-8").decode(buffer8.slice(old_path_ptr, old_path_ptr + old_path_len));
					const new_path = new TextDecoder("utf-8").decode(buffer8.slice(new_path_ptr, new_path_ptr + new_path_len));
					let { ret, inode_obj } = self.fds[fd].path_unlink(old_path);
					if (inode_obj == null) return ret;
					ret = self.fds[new_fd].path_link(new_path, inode_obj, true);
					if (ret != 0) {
						if (self.fds[fd].path_link(old_path, inode_obj, true) != 0) throw "path_link should always return success when relinking an inode back to the original place";
					}
					return ret;
				} else return 8;
			},
			path_symlink(old_path_ptr, old_path_len, fd, new_path_ptr, new_path_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					new TextDecoder("utf-8").decode(buffer8.slice(old_path_ptr, old_path_ptr + old_path_len));
					new TextDecoder("utf-8").decode(buffer8.slice(new_path_ptr, new_path_ptr + new_path_len));
					return 58;
				} else return 8;
			},
			path_unlink_file(fd, path_ptr, path_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer);
				if (self.fds[fd] != void 0) {
					const path = new TextDecoder("utf-8").decode(buffer8.slice(path_ptr, path_ptr + path_len));
					return self.fds[fd].path_unlink_file(path);
				} else return 8;
			},
			poll_oneoff(in_ptr, out_ptr, nsubscriptions) {
				if (nsubscriptions === 0) return 28;
				if (nsubscriptions > 1) {
					debug.log("poll_oneoff: only a single subscription is supported");
					return 58;
				}
				const buffer = new DataView(self.inst.exports.memory.buffer);
				const s = Subscription.read_bytes(buffer, in_ptr);
				const eventtype = s.eventtype;
				const clockid = s.clockid;
				const timeout = s.timeout;
				if (eventtype !== 0) {
					debug.log("poll_oneoff: only clock subscriptions are supported");
					return 58;
				}
				let getNow = void 0;
				if (clockid === 1) getNow = () => BigInt(Math.round(performance.now() * 1e6));
				else if (clockid === 0) getNow = () => BigInt((/* @__PURE__ */ new Date()).getTime()) * 1000000n;
				else return 28;
				const endTime = (s.flags & 1) !== 0 ? timeout : getNow() + timeout;
				while (endTime > getNow());
				new Event(s.userdata, 0, eventtype).write_bytes(buffer, out_ptr);
				return 0;
			},
			proc_exit(exit_code) {
				throw new WASIProcExit(exit_code);
			},
			proc_raise(sig) {
				throw "raised signal " + sig;
			},
			sched_yield() {},
			random_get(buf, buf_len) {
				const buffer8 = new Uint8Array(self.inst.exports.memory.buffer).subarray(buf, buf + buf_len);
				if ("crypto" in globalThis && (typeof SharedArrayBuffer === "undefined" || !(self.inst.exports.memory.buffer instanceof SharedArrayBuffer))) for (let i = 0; i < buf_len; i += 65536) crypto.getRandomValues(buffer8.subarray(i, i + 65536));
				else for (let i = 0; i < buf_len; i++) buffer8[i] = Math.random() * 256 | 0;
			},
			sock_recv(fd, ri_data, ri_flags) {
				throw "sockets not supported";
			},
			sock_send(fd, si_data, si_flags) {
				throw "sockets not supported";
			},
			sock_shutdown(fd, how) {
				throw "sockets not supported";
			},
			sock_accept(fd, flags) {
				throw "sockets not supported";
			}
		};
	}
};
//#endregion
//#region ../swig-flow-refactor/node_modules/@bjorn3/browser_wasi_shim/dist/fd.js
var Inode = class Inode {
	static issue_ino() {
		return Inode.next_ino++;
	}
	static root_ino() {
		return 0n;
	}
	constructor() {
		this.ino = Inode.issue_ino();
	}
};
Inode.next_ino = 1n;
//#endregion
//#region src/balanced-calling-profile.ts
const decoder$1 = new TextDecoder();
const encoder$1 = new TextEncoder();
function longestExactRun(query, germline) {
	let current = 0;
	let longest = 0;
	for (let index = 0; index < Math.min(query.length, germline.length); index += 1) if (query[index] !== "-" && query[index] === germline[index]) {
		current += 1;
		longest = Math.max(longest, current);
	} else current = 0;
	return longest;
}
const complement = {
	A: "T",
	C: "G",
	G: "C",
	T: "A",
	U: "A",
	R: "Y",
	Y: "R",
	M: "K",
	K: "M",
	S: "S",
	W: "W",
	B: "V",
	V: "B",
	D: "H",
	H: "D",
	N: "N",
	"-": "-"
};
function reverseComplement(sequence) {
	let result = "";
	for (let index = sequence.length - 1; index >= 0; index -= 1) result += complement[sequence[index].toUpperCase()] ?? "N";
	return result;
}
const aminoAcids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const baseIndex = {
	T: 0,
	C: 1,
	A: 2,
	G: 3
};
function translateAlignment(sequence, frame) {
	let protein = "";
	for (let index = frame; index + 2 < sequence.length; index += 3) {
		const codon = sequence.slice(index, index + 3).toUpperCase();
		if (codon.includes("-")) {
			protein += "-";
			continue;
		}
		const a = baseIndex[codon[0]];
		const b = baseIndex[codon[1]];
		const c = baseIndex[codon[2]];
		protein += a === void 0 || b === void 0 || c === void 0 ? "X" : aminoAcids[a * 16 + b * 4 + c];
	}
	return protein;
}
function indexColumns(header) {
	return new Map(header.split("	").map((name, index) => [name, index]));
}
function field(fields, columns, name) {
	const index = columns.get(name);
	return index === void 0 ? "" : fields[index] ?? "";
}
function setField(fields, columns, name, value = "") {
	const index = columns.get(name);
	if (index !== void 0) fields[index] = value;
}
function rebuildCompositeAlignment(fields, columns) {
	const sequence = field(fields, columns, "sequence");
	const oriented = field(fields, columns, "rev_comp") === "T" ? reverseComplement(sequence) : sequence;
	const vEnd = Number(field(fields, columns, "v_sequence_end"));
	const jStart = Number(field(fields, columns, "j_sequence_start"));
	const gap = Number.isFinite(vEnd) && Number.isFinite(jStart) && jStart > vEnd ? oriented.slice(vEnd, jStart - 1) : "";
	const query = field(fields, columns, "v_sequence_alignment") + gap + field(fields, columns, "j_sequence_alignment");
	const germline = field(fields, columns, "v_germline_alignment") + "N".repeat(gap.length) + field(fields, columns, "j_germline_alignment");
	setField(fields, columns, "sequence_alignment", query);
	setField(fields, columns, "germline_alignment", germline);
	setField(fields, columns, "np1", gap);
	setField(fields, columns, "np1_length", String(Math.max(0, jStart - 1 - vEnd)));
	setField(fields, columns, "np2");
	setField(fields, columns, "np2_length");
	const vStart = Number(field(fields, columns, "v_sequence_start"));
	const cdr3Start = Number(field(fields, columns, "cdr3_start"));
	let cysColumn = -1;
	if (Number.isFinite(vStart) && Number.isFinite(cdr3Start)) {
		const cysPosition = cdr3Start - 4;
		let queryPosition = vStart - 1;
		for (let column = 0; column < query.length; column += 1) {
			if (query[column] === "-") continue;
			if (queryPosition === cysPosition) {
				cysColumn = column;
				break;
			}
			queryPosition += 1;
		}
	}
	if (cysColumn >= 0) {
		const frame = cysColumn % 3;
		setField(fields, columns, "sequence_alignment_aa", translateAlignment(query, frame));
		setField(fields, columns, "germline_alignment_aa", translateAlignment(germline, frame));
	} else {
		setField(fields, columns, "sequence_alignment_aa");
		setField(fields, columns, "germline_alignment_aa");
	}
}
/**
* Apply the calibrated IgBLAST-balanced evidence rule to an AIRR batch.
*
* Only D hits whose strongest alignment has exactly a five-nucleotide exact
* run are candidates. Such a hit is removed when
* `j_sequence_start - v_sequence_end <= 11`; all longer exact runs and all
* five-base hits in longer V-J spans are retained. D-dependent AIRR fields
* and the stitched V-J alignment are rebuilt so the row stays internally
* consistent.
*/
function applyBalancedDFilter(header, body) {
	const columns = indexColumns(header);
	if ([
		"sequence_id",
		"d_call",
		"d_sequence_alignment",
		"d_germline_alignment",
		"v_sequence_end",
		"j_sequence_start"
	].some((name) => !columns.has(name))) throw new Error("The AIRR table is missing fields required by the IgBLAST-balanced D rule.");
	const text = decoder$1.decode(body);
	const hadTrailingNewline = text.endsWith("\n");
	const lines = text.split("\n");
	if (hadTrailingNewline) lines.pop();
	const suppressedSequenceIds = /* @__PURE__ */ new Set();
	const dFields = [
		"d_call",
		"d_score",
		"d_identity",
		"d_cigar",
		"d_sequence_start",
		"d_sequence_end",
		"d_germline_start",
		"d_germline_end",
		"d_sequence_alignment",
		"d_germline_alignment",
		"d_frame",
		"d_alternatives"
	];
	const transformed = lines.map((line) => {
		if (!line) return line;
		const fields = line.replace(/\r$/, "").split("	");
		if (!field(fields, columns, "d_call")) return line.replace(/\r$/, "");
		const exact = longestExactRun(field(fields, columns, "d_sequence_alignment"), field(fields, columns, "d_germline_alignment"));
		const vjSpan = Math.max(0, Number(field(fields, columns, "j_sequence_start")) - Number(field(fields, columns, "v_sequence_end")));
		if (exact !== 5 || vjSpan > 11) return line.replace(/\r$/, "");
		suppressedSequenceIds.add(field(fields, columns, "sequence_id"));
		dFields.forEach((name) => setField(fields, columns, name));
		rebuildCompositeAlignment(fields, columns);
		return fields.join("	");
	});
	if (!suppressedSequenceIds.size) return {
		body,
		suppressedSequenceIds
	};
	return {
		body: encoder$1.encode(transformed.join("\n") + (hadTrailingNewline ? "\n" : "")),
		suppressedSequenceIds
	};
}
function reconcileBalancedDoubleD(header, body, suppressedSequenceIds) {
	if (!suppressedSequenceIds.size) return body;
	const columns = indexColumns(header);
	const idIndex = columns.get("sequence_id");
	const standardDIndex = columns.get("standard_d_call");
	if (idIndex === void 0 || standardDIndex === void 0) return body;
	const text = decoder$1.decode(body);
	const hadTrailingNewline = text.endsWith("\n");
	const lines = text.split("\n");
	if (hadTrailingNewline) lines.pop();
	let changed = false;
	const transformed = lines.map((line) => {
		const fields = line.replace(/\r$/, "").split("	");
		if (!suppressedSequenceIds.has(fields[idIndex] ?? "")) return line.replace(/\r$/, "");
		fields[standardDIndex] = "";
		changed = true;
		return fields.join("	");
	});
	return changed ? encoder$1.encode(transformed.join("\n") + (hadTrailingNewline ? "\n" : "")) : body;
}
//#endregion
//#region cli-src/swig-worker.mjs
const encoder = new TextEncoder();
const decoder = new TextDecoder();
let runtime = null;
let callingProfile = "truth_optimized";
function put(bytes) {
	const pointer = runtime.swig_alloc(bytes.byteLength);
	if (!pointer && bytes.byteLength) throw new Error("SwiftIG ran out of WebAssembly memory.");
	new Uint8Array(runtime.memory.buffer, pointer, bytes.byteLength).set(bytes);
	return [pointer, bytes.byteLength];
}
function read(pointer, length) {
	return decoder.decode(new Uint8Array(runtime.memory.buffer, pointer, length));
}
function errorText() {
	return read(runtime.swig_error_ptr(), runtime.swig_error_len()) || "SwiftIG could not complete the annotation.";
}
async function initialize(message) {
	const wasi = new WASI([], [], []);
	const module = await WebAssembly.compile(readFileSync(message.wasmPath));
	const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
	wasi.initialize(instance);
	runtime = instance.exports;
	const strategy = message.assignerStrategy === "standard" ? 0 : message.assignerStrategy === "aer" ? 2 : 1;
	if (runtime.swig_set_assigner_strategy(strategy) !== 0) throw new Error("SwiftIG rejected the assignment strategy.");
	const profile = message.callingProfile === "truth_optimized" ? 0 : 1;
	if (runtime.swig_set_calling_profile(profile) !== 0) throw new Error("SwiftIG rejected the calling profile.");
	callingProfile = message.callingProfile;
	if (message.hasTuning) {
		if (typeof runtime.swig_set_tuning_options !== "function") throw new Error("This SwiftIG build does not expose the requested D/J compatibility controls.");
		if (runtime.swig_set_tuning_options(message.tuningDMatch, message.tuningDMismatch, message.tuningDGapOpen, message.tuningDGapExtend, message.tuningTopD, message.tuningMinDMatch, message.tuningJMatch, message.tuningJMismatch, message.tuningJGapOpen, message.tuningJGapExtend, message.tuningTopJ, message.tuningMinJLength) !== 0) throw new Error("SwiftIG rejected the requested D/J compatibility controls.");
	}
	const allocations = [
		message.referenceV,
		message.referenceD,
		message.referenceJ,
		message.referenceC
	].map((value) => put(encoder.encode(value || "")));
	try {
		const genes = runtime.swig_init_database(...allocations.flat());
		if (genes < 0) throw new Error(errorText());
		return { genes };
	} finally {
		allocations.forEach(([pointer]) => runtime.swig_free(pointer));
	}
}
function annotate(message) {
	const [pointer, size] = put(encoder.encode(message.text));
	let count;
	try {
		count = message.doubleDMode !== "off" ? runtime.swig_annotate_double_d(pointer, size, message.format, Math.round(message.minimumIdentity * 1e3), message.strand, message.doubleDMode === "all" ? 1 : 2, Math.round(message.doubleDMinimumVjSpan), Math.round(message.doubleDSeedLength), Math.round(message.doubleDPseudoTrim), Math.round(message.doubleDMaximumPseudoMismatches), Math.round(message.doubleDMinimumScoreGain)) : runtime.swig_annotate(pointer, size, message.format, Math.round(message.minimumIdentity * 1e3), message.strand);
	} finally {
		runtime.swig_free(pointer);
	}
	if (count < 0) throw new Error(errorText());
	if (count !== message.count) throw new Error(`SwiftIG returned ${count} rows for a ${message.count}-record batch.`);
	const result = new Uint8Array(runtime.memory.buffer, runtime.swig_result_ptr(), runtime.swig_result_len());
	const newline = result.indexOf(10);
	if (newline < 0) throw new Error("SwiftIG returned an invalid AIRR table.");
	const header = decoder.decode(result.subarray(0, newline)).replace(/\r$/, "");
	let body = result.subarray(newline + 1);
	const balanced = callingProfile === "igblast_balanced" ? applyBalancedDFilter(header, body) : null;
	if (balanced) body = balanced.body;
	const response = {
		header,
		body: decoder.decode(body),
		count,
		doubleDHeader: "",
		doubleDBody: ""
	};
	if (message.doubleDMode !== "off") {
		const dd = new Uint8Array(runtime.memory.buffer, runtime.swig_double_d_result_ptr(), runtime.swig_double_d_result_len());
		const ddNewline = dd.indexOf(10);
		if (ddNewline < 0) throw new Error("SwiftIG returned invalid double-D evidence.");
		response.doubleDHeader = decoder.decode(dd.subarray(0, ddNewline)).replace(/\r$/, "");
		let ddBody = dd.slice(ddNewline + 1);
		if (balanced) ddBody = reconcileBalancedDoubleD(response.doubleDHeader, ddBody, balanced.suppressedSequenceIds);
		response.doubleDBody = decoder.decode(ddBody);
	}
	return response;
}
const send = (message) => parentPort ? parentPort.postMessage(message) : globalThis.postMessage(message);
const receive = async (message) => {
	try {
		const result = message.type === "init" ? await initialize(message) : annotate(message);
		send({
			id: message.id,
			...result
		});
	} catch (error) {
		send({
			id: message.id,
			error: error instanceof Error ? error.message : String(error)
		});
	}
};
if (parentPort) parentPort.on("message", receive);
else globalThis.addEventListener("message", (event) => void receive(event.data));
//#endregion
